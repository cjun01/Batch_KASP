#!/usr/bin/env python3
"""Check whether marker CSV coordinates and ref alleles match a reference FASTA."""

import argparse
import csv
from collections import defaultdict

try:
    from Bio import SeqIO
except ImportError:  # pragma: no cover - fallback for lightweight environments
    SeqIO = None


def normalize_chrom_name(name):
    value = str(name).strip()
    lowered = value.lower()
    if lowered.startswith("chr"):
        return lowered[3:]
    return lowered


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Validate a marker CSV against a reference FASTA. "
            "The marker CSV must contain Chr, position, and ref columns."
        )
    )
    parser.add_argument(
        "-g",
        "--genome_path",
        required=True,
        help="Path to the reference genome FASTA file.",
    )
    parser.add_argument(
        "-m",
        "--marker_csv",
        required=True,
        help="Path to the marker CSV file.",
    )
    parser.add_argument(
        "-o",
        "--output_csv",
        help=(
            "Optional path for a per-marker validation report. "
            "If omitted, only the summary is printed."
        ),
    )
    return parser.parse_args()


def load_markers(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        missing = {"Chr", "position", "ref"} - set(reader.fieldnames or [])
        if missing:
            missing_columns = ", ".join(sorted(missing))
            raise ValueError(
                "Marker CSV must contain Chr, position, and ref columns. "
                f"Missing: {missing_columns}"
            )

        rows = []
        grouped = defaultdict(list)
        for index, row in enumerate(reader, start=1):
            chrom = row["Chr"]
            position = int(row["position"])
            ref = row["ref"].upper()
            alt = row.get("alt", "")
            normalized_chrom = normalize_chrom_name(chrom)
            record = {
                "row_index": index,
                "Chr": chrom,
                "position": position,
                "ref": ref,
                "alt": alt.upper() if alt else "",
                "status": "missing_chromosome_in_fasta",
                "fasta_base": "",
            }
            rows.append(record)
            grouped[normalized_chrom].append(record)

    return rows, grouped


def check_markers(rows_by_chrom, genome_path):
    found_chromosomes = set()

    if SeqIO is not None:
        with open(genome_path, "r") as handle:
            record_iter = (
                (record.id, str(record.seq).upper())
                for record in SeqIO.parse(handle, "fasta")
            )
            process_fasta_records(record_iter, rows_by_chrom, found_chromosomes)
    else:
        with open(genome_path, "r") as handle:
            process_fasta_records(simple_fasta_iter(handle), rows_by_chrom, found_chromosomes)

    for chrom_key, rows in rows_by_chrom.items():
        if chrom_key in found_chromosomes:
            continue
        for row in rows:
            row["status"] = "missing_chromosome_in_fasta"
            row["fasta_base"] = ""


def process_fasta_records(record_iter, rows_by_chrom, found_chromosomes):
    for record_id, sequence in record_iter:
        chrom_key = normalize_chrom_name(record_id)
        if chrom_key not in rows_by_chrom:
            continue

        found_chromosomes.add(chrom_key)
        sequence_length = len(sequence)

        for row in rows_by_chrom[chrom_key]:
            pos = row["position"]
            if pos < 1 or pos > sequence_length:
                row["status"] = "position_out_of_range"
                row["fasta_base"] = ""
                continue

            fasta_base = sequence[pos - 1]
            row["fasta_base"] = fasta_base
            row["status"] = "match" if fasta_base == row["ref"] else "ref_mismatch"


def simple_fasta_iter(handle):
    current_id = None
    sequence_parts = []

    for line in handle:
        if line.startswith(">"):
            if current_id is not None:
                yield current_id, "".join(sequence_parts).upper()
            current_id = line[1:].strip().split()[0]
            sequence_parts = []
        else:
            sequence_parts.append(line.strip())

    if current_id is not None:
        yield current_id, "".join(sequence_parts).upper()


def write_report(rows, path):
    fieldnames = [
        "row_index",
        "Chr",
        "position",
        "ref",
        "alt",
        "fasta_base",
        "status",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    rows, rows_by_chrom = load_markers(args.marker_csv)
    check_markers(rows_by_chrom, args.genome_path)

    counts = defaultdict(int)
    for row in rows:
        counts[row["status"]] += 1

    print(f"Total markers checked: {len(rows)}")
    print(f"Matches: {counts['match']}")
    print(f"Reference mismatches: {counts['ref_mismatch']}")
    print(f"Missing chromosomes: {counts['missing_chromosome_in_fasta']}")
    print(f"Positions out of range: {counts['position_out_of_range']}")

    if args.output_csv:
        write_report(rows, args.output_csv)
        print(f"Detailed report written to {args.output_csv}")


if __name__ == "__main__":
    main()
