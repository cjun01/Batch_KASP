#!/usr/bin/env python3
"""Convert a VCF into the marker CSV expected by KASP_design.py."""

import argparse
import csv
import gzip
import sys


VALID_BASES = {"A", "C", "G", "T"}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Convert a VCF/VCF.gz into the marker CSV expected by "
            "KASP_design.py (Chr,position,ref,alt)."
        )
    )
    parser.add_argument(
        "-i",
        "--input_vcf",
        required=True,
        help="Path to the input VCF or VCF.gz file.",
    )
    parser.add_argument(
        "-o",
        "--output_csv",
        required=True,
        help="Path to the output CSV file.",
    )
    parser.add_argument(
        "-s",
        "--sample",
        help=(
            "Optional sample name to filter variants by genotype. "
            "When provided, --genotype defaults to 'alt'."
        ),
    )
    parser.add_argument(
        "--genotype",
        choices=("alt", "nonref", "any"),
        help=(
            "Genotype filter to apply when --sample is provided: "
            "'alt' keeps only 1/1 calls, 'nonref' keeps 0/1 and 1/1, "
            "and 'any' keeps any non-missing call."
        ),
    )
    parser.add_argument(
        "--exclude-contig-substring",
        default=None,
        help="Skip contigs whose names contain this substring, e.g. 'unitig'.",
    )
    return parser.parse_args()


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", newline="")


def normalize_genotype(raw_gt):
    gt = raw_gt.replace("|", "/")
    if gt == "." or "." in gt:
        return "missing"

    alleles = gt.split("/")
    if len(alleles) != 2:
        return "other"
    if any(allele not in {"0", "1"} for allele in alleles):
        return "other"
    if alleles == ["0", "0"]:
        return "ref"
    if alleles == ["1", "1"]:
        return "alt"
    if set(alleles) == {"0", "1"}:
        return "het"
    return "other"


def extract_genotype(format_field, sample_field):
    format_keys = format_field.split(":") if format_field else []
    sample_values = sample_field.split(":") if sample_field else []

    if "GT" in format_keys:
        gt_index = format_keys.index("GT")
        if gt_index >= len(sample_values):
            return "missing"
        raw_gt = sample_values[gt_index]
    elif sample_values:
        raw_gt = sample_values[0]
    else:
        return "missing"

    return normalize_genotype(raw_gt)


def genotype_allowed(state, genotype_mode):
    if genotype_mode == "alt":
        return state == "alt"
    if genotype_mode == "nonref":
        return state in {"het", "alt"}
    if genotype_mode == "any":
        return state in {"ref", "het", "alt"}
    return True


def is_simple_snp(ref, alt):
    if "," in alt:
        return False
    if len(ref) != 1 or len(alt) != 1:
        return False
    return ref.upper() in VALID_BASES and alt.upper() in VALID_BASES


def main():
    args = parse_args()
    genotype_mode = args.genotype or ("alt" if args.sample else None)

    stats = {
        "records_seen": 0,
        "written": 0,
        "non_simple_snp": 0,
        "genotype_filtered": 0,
        "excluded_contig": 0,
        "duplicate": 0,
    }
    seen = set()
    header = None
    sample_index = None

    with open_text(args.input_vcf) as infile, open(args.output_csv, "w", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=["Chr", "position", "ref", "alt"])
        writer.writeheader()

        for raw_line in infile:
            if raw_line.startswith("##"):
                continue

            line = raw_line.rstrip("\n")
            if not line:
                continue

            if line.startswith("#CHROM"):
                header = line.split("\t")
                if args.sample:
                    samples = header[9:]
                    if args.sample not in samples:
                        available = ", ".join(samples[:10])
                        suffix = "" if len(samples) <= 10 else ", ..."
                        raise SystemExit(
                            f"Sample '{args.sample}' was not found in the VCF header. "
                            f"Available samples: {available}{suffix}"
                        )
                    sample_index = header.index(args.sample)
                continue

            if header is None:
                raise SystemExit("VCF header line starting with '#CHROM' was not found.")

            fields = line.split("\t")
            if len(fields) < 8:
                raise SystemExit(f"Malformed VCF row with fewer than 8 columns: {line}")

            stats["records_seen"] += 1
            chrom, pos, _variant_id, ref, alt = fields[:5]

            if args.exclude_contig_substring and args.exclude_contig_substring in chrom:
                stats["excluded_contig"] += 1
                continue

            ref = ref.upper()
            alt = alt.upper()
            if not is_simple_snp(ref, alt):
                stats["non_simple_snp"] += 1
                continue

            if sample_index is not None:
                if len(fields) <= sample_index:
                    raise SystemExit(
                        f"VCF row is missing sample column '{args.sample}': {line}"
                    )
                format_field = fields[8] if len(fields) > 8 else ""
                genotype_state = extract_genotype(format_field, fields[sample_index])
                if not genotype_allowed(genotype_state, genotype_mode):
                    stats["genotype_filtered"] += 1
                    continue

            key = (chrom, pos, ref, alt)
            if key in seen:
                stats["duplicate"] += 1
                continue
            seen.add(key)

            writer.writerow(
                {
                    "Chr": chrom,
                    "position": pos,
                    "ref": ref,
                    "alt": alt,
                }
            )
            stats["written"] += 1

    print(f"Wrote {stats['written']} markers to {args.output_csv}")
    print(f"Records processed: {stats['records_seen']}")
    if args.sample:
        print(f"Sample filter: {args.sample} ({genotype_mode})")
    print(f"Skipped non-simple SNPs: {stats['non_simple_snp']}")
    print(f"Skipped by genotype filter: {stats['genotype_filtered']}")
    print(f"Skipped by contig filter: {stats['excluded_contig']}")
    print(f"Skipped duplicates: {stats['duplicate']}")

    if stats["written"] == 0:
        print(
            "Warning: no markers were written. Check the sample name, genotype filter, "
            "or contig filter.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
