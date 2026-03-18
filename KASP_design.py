import argparse
import gzip
from bisect import bisect_left, bisect_right
from collections import defaultdict

from Bio import SeqIO
import pandas as pd

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover - fallback when tqdm is not installed
    def tqdm(iterable, **_kwargs):
        return iterable


FAM = "GAAGGTGACCAAGTTCATGCT"
VIC = "GAAGGTCGGAGTCAACGGATT"
VALID_BASES = {"A", "C", "G", "T"}


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def normalize_chrom_name(name):
    value = str(name).strip()
    lowered = value.lower()
    if lowered.startswith("chr"):
        return lowered[3:]
    return lowered


def add_variant_span(background_positions, chrom, start, span=1):
    chrom_key = normalize_chrom_name(chrom)
    for position in range(start, start + max(1, span)):
        background_positions[chrom_key].add(position)


def normalize_marker_columns(df):
    column_map = {column.lower(): column for column in df.columns}
    required = {"chr", "position", "ref", "alt"}
    missing = required - set(column_map)
    if missing:
        missing_columns = ", ".join(sorted(missing))
        raise ValueError(
            "Marker CSV must contain the columns Chr, position, ref, and alt. "
            f"Missing: {missing_columns}"
        )

    normalized = pd.DataFrame(
        {
            "Chr": df[column_map["chr"]].astype(str),
            "position": pd.to_numeric(df[column_map["position"]], errors="raise").astype(int),
            "ref": df[column_map["ref"]].astype(str).str.upper(),
            "alt": df[column_map["alt"]].astype(str).str.upper(),
        }
    )

    valid_mask = (
        normalized["ref"].str.len().eq(1)
        & normalized["alt"].str.len().eq(1)
        & normalized["ref"].isin(VALID_BASES)
        & normalized["alt"].isin(VALID_BASES)
        & normalized["ref"].ne(normalized["alt"])
    )

    skipped = int((~valid_mask).sum())
    if skipped:
        print(f"Skipping {skipped} non-SNP or invalid markers from the marker CSV.")

    normalized = normalized.loc[valid_mask].drop_duplicates().sort_values(["Chr", "position"])
    return normalized.reset_index(drop=True)


def load_marker_table(path):
    return normalize_marker_columns(pd.read_csv(path))


def load_background_from_csv(path, background_positions):
    df = pd.read_csv(path)
    column_map = {column.lower(): column for column in df.columns}
    required = {"chr", "position"}
    missing = required - set(column_map)
    if missing:
        missing_columns = ", ".join(sorted(missing))
        raise ValueError(
            "Background CSV must contain at least the columns Chr and position. "
            f"Missing: {missing_columns}"
        )

    ref_column = column_map.get("ref")
    positions = pd.to_numeric(df[column_map["position"]], errors="raise").astype(int)
    chromosomes = df[column_map["chr"]].astype(str)

    for index, chrom in chromosomes.items():
        start = int(positions.loc[index])
        if ref_column is not None:
            ref_value = str(df.at[index, ref_column]).upper()
            span = len(ref_value) if ref_value not in {"", "NAN"} else 1
        else:
            span = 1
        add_variant_span(background_positions, chrom, start, span)


def load_background_from_vcf(path, background_positions):
    with open_text(path) as handle:
        for raw_line in handle:
            if not raw_line or raw_line.startswith("#"):
                continue

            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue

            chrom = fields[0]
            start = int(fields[1])
            ref = fields[3].upper()
            span = len(ref) if ref else 1
            add_variant_span(background_positions, chrom, start, span)


def build_background_positions(targets, background_csv=None, background_vcf=None):
    background_positions = defaultdict(set)

    # Always include the target set so nearby target markers are also avoided.
    for row in targets.itertuples(index=False):
        add_variant_span(background_positions, row.Chr, int(row.position), len(row.ref))

    if background_csv:
        load_background_from_csv(background_csv, background_positions)
    if background_vcf:
        load_background_from_vcf(background_vcf, background_positions)

    return {
        chrom: sorted(positions)
        for chrom, positions in background_positions.items()
    }


class MarkerDesign:
    def __init__(self, background_positions=None):
        self.background_positions = background_positions or {}

    def get_complement(self, sequence):
        complements = {"A": "T", "T": "A", "C": "G", "G": "C"}
        return "".join(complements[base] for base in sequence)

    def reverse_complement(self, seq):
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        return "".join(complement[base] for base in reversed(seq))

    def find_TM(self, seq):
        num_g = seq.count("G")
        num_c = seq.count("C")
        num_a = seq.count("A")
        num_t = seq.count("T")
        return 64.9 + 41 * (num_g + num_c - 16.4) / (num_g + num_c + num_a + num_t)

    def primer_dimer_check(self, primer1, primer2, threshold=5):
        end_primer1 = primer1[-threshold:]
        end_primer2 = self.reverse_complement(primer2[-threshold:])
        return end_primer1 == end_primer2

    def find_hairpins(self, seq, min_loop=3, min_stem=4):
        for index in range(len(seq)):
            for end in range(index + min_stem * 2 + min_loop, len(seq) + 1):
                stem1 = seq[index:index + min_stem]
                stem2 = seq[end - min_stem:end]
                complement_stem2 = self.get_complement(stem2)
                if stem1 == complement_stem2:
                    return True
        return False

    def calculate_gc_content(self, sequence):
        sequence = sequence.upper()
        gc_count = sequence.count("G") + sequence.count("C")
        return (gc_count / len(sequence)) * 100

    def find_repeats(self, sequence):
        thresholds = {1: 4, 2: 4, 3: 4}
        sequence_length = len(sequence)

        for motif_length, min_repeats in thresholds.items():
            for index in range(sequence_length - motif_length + 1):
                repeat = sequence[index:index + motif_length]
                max_repeats = 0
                for multiplier in range(1, (sequence_length - index) // motif_length + 1):
                    if sequence[index:index + motif_length * multiplier] == repeat * multiplier:
                        max_repeats = multiplier
                    else:
                        break
                if max_repeats >= min_repeats:
                    return (
                        "Yes",
                        f"{max_repeats} repetitions of a {motif_length}-nucleotide repeat: '{repeat}'",
                    )
        return "No"

    def check(self, sequence, position_zero_based, ref):
        if sequence is None:
            return False
        if position_zero_based < 0 or position_zero_based >= len(sequence):
            return False
        return sequence[position_zero_based] == ref

    def region_has_background_variant(self, chrom, start, end, ignore_positions=None):
        chrom_key = normalize_chrom_name(chrom)
        positions = self.background_positions.get(chrom_key)
        if not positions:
            return False

        left = bisect_left(positions, start)
        right = bisect_right(positions, end)
        if left == right:
            return False

        if not ignore_positions:
            return True

        ignored = set(ignore_positions)
        for position in positions[left:right]:
            if position not in ignored:
                return True
        return False

    def find_reverse_primers(
        self,
        chrom,
        sequence,
        flank_start_one_based,
        target_tm,
        target_position_one_based,
        tm_tolerance=2,
        min_length=20,
        max_length=30,
        max_candidates=5,
    ):
        potential_primers = []
        summary = {
            "insufficient_flank": False,
            "tm_mismatch": 0,
            "repeat_fail": 0,
            "gc_fail": 0,
            "background_fail": 0,
        }
        if len(sequence) < min_length:
            summary["insufficient_flank"] = True
            return potential_primers, summary

        for start_offset in range(0, len(sequence) - min_length + 1):
            for primer_length in range(min_length, max_length + 1):
                end_offset = start_offset + primer_length
                if end_offset > len(sequence):
                    break

                binding_sequence = sequence[start_offset:end_offset]
                primer = self.reverse_complement(binding_sequence)
                reverse_tm = self.find_TM(primer)

                if not (target_tm - tm_tolerance <= reverse_tm <= target_tm + tm_tolerance):
                    summary["tm_mismatch"] += 1
                    continue
                if self.find_repeats(primer) != "No":
                    summary["repeat_fail"] += 1
                    continue
                if self.calculate_gc_content(primer) < 40:
                    summary["gc_fail"] += 1
                    continue

                primer_start = flank_start_one_based + start_offset
                primer_end = flank_start_one_based + end_offset - 1
                if self.region_has_background_variant(
                    chrom,
                    primer_start,
                    primer_end,
                    ignore_positions={target_position_one_based},
                ):
                    summary["background_fail"] += 1
                    continue

                potential_primers.append((primer, reverse_tm, primer_start, primer_end))
                if len(potential_primers) >= max_candidates:
                    return potential_primers, summary

        return potential_primers, summary

    def design_kasp(self, sequences, chrom, position, ref, alt, reverse_window=200):
        sequence = sequences.get(normalize_chrom_name(chrom))
        if sequence is None:
            return None, "missing_chromosome_in_fasta", False

        target_position_one_based = int(position)
        target_position_zero_based = target_position_one_based - 1
        if target_position_zero_based < 0 or target_position_zero_based >= len(sequence):
            return None, "position_out_of_range", False

        ref = ref.upper()
        alt = alt.upper()
        alignment_check = self.check(sequence, target_position_zero_based, ref)

        right_flank_start_zero_based = target_position_zero_based + 1
        right_flank = sequence[
            right_flank_start_zero_based:right_flank_start_zero_based + reverse_window
        ]

        found_tm_match = False
        forward_blocked_by_variant = False
        reverse_failure_reason = "no_reverse_primer_candidate"

        for forward_core_length in range(18, 36):
            forward_start_zero_based = target_position_zero_based - (forward_core_length - 1)
            if forward_start_zero_based < 0:
                continue

            forward_end_zero_based = target_position_zero_based
            forward_template = sequence[
                forward_start_zero_based:forward_end_zero_based + 1
            ]
            forward_ref_core = forward_template[:-1] + ref
            forward_alt_core = forward_template[:-1] + alt
            tm_forward = self.find_TM(forward_ref_core)

            if not 57 <= tm_forward <= 65:
                continue
            found_tm_match = True
            if self.region_has_background_variant(
                chrom,
                forward_start_zero_based + 1,
                forward_end_zero_based + 1,
                ignore_positions={target_position_one_based},
            ):
                forward_blocked_by_variant = True
                continue

            a1_primer = FAM + forward_ref_core
            a2_primer = VIC + forward_alt_core
            reverse_candidates, reverse_summary = self.find_reverse_primers(
                chrom,
                right_flank,
                right_flank_start_zero_based + 1,
                tm_forward,
                target_position_one_based,
            )
            if not reverse_candidates:
                if reverse_summary["insufficient_flank"]:
                    reverse_failure_reason = "insufficient_downstream_flank"
                elif reverse_summary["background_fail"] > 0:
                    reverse_failure_reason = "reverse_primer_overlaps_variant"
                elif reverse_summary["gc_fail"] > 0:
                    reverse_failure_reason = "reverse_primer_gc_filter"
                elif reverse_summary["repeat_fail"] > 0:
                    reverse_failure_reason = "reverse_primer_repeat_filter"
                else:
                    reverse_failure_reason = "reverse_primer_tm_filter"

            best_fallback = None
            for reverse_primer, reverse_tm, reverse_start, reverse_end in reverse_candidates:
                product_size = reverse_end - (forward_start_zero_based + 1) + 1
                has_hairpin = (
                    self.find_hairpins(forward_ref_core)
                    or self.find_hairpins(forward_alt_core)
                    or self.find_hairpins(reverse_primer)
                )
                has_dimer = (
                    self.primer_dimer_check(a1_primer, reverse_primer)
                    or self.primer_dimer_check(a2_primer, reverse_primer)
                )

                result = {
                    "A1_Primer": a1_primer,
                    "A2_Primer": a2_primer,
                    "Reverse_Primer": reverse_primer,
                    "Tm": reverse_tm,
                    "Product size": product_size,
                    "Hairpin": "Yes" if has_hairpin else "No",
                    "Primer dimer": "Yes" if has_dimer else "No",
                    "GC content for forward primer": self.calculate_gc_content(a1_primer),
                    "Repeats within forward primer": self.find_repeats(a1_primer),
                    "GC content for reverse primer": self.calculate_gc_content(reverse_primer),
                    "Repeats within reverse primer": self.find_repeats(reverse_primer),
                    "SNP alignment Check": alignment_check,
                }

                if not has_hairpin and not has_dimer:
                    return result, "designed", alignment_check
                if best_fallback is None:
                    best_fallback = result

            if best_fallback is not None:
                return best_fallback, "fallback_hairpin_or_dimer", alignment_check

        if not found_tm_match:
            return None, "no_forward_tm_match", alignment_check
        if forward_blocked_by_variant:
            return None, "forward_primer_overlaps_variant", alignment_check
        return None, reverse_failure_reason, alignment_check


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate primers for KASP marker design.")
    parser.add_argument(
        "-g",
        "--genome_path",
        required=True,
        help="Path to the genome FASTA file.",
    )
    parser.add_argument(
        "-m",
        "--marker_csv",
        required=True,
        help="Path to the marker CSV containing Chr, position, ref, and alt.",
    )
    parser.add_argument(
        "-o",
        "--output_csv",
        required=True,
        help="Path to save the output CSV file with primer designs.",
    )
    parser.add_argument(
        "--background_vcf",
        help=(
            "Optional VCF or VCF.gz of background variants to avoid inside primer "
            "binding regions."
        ),
    )
    parser.add_argument(
        "--background_csv",
        help=(
            "Optional background CSV of variants to avoid inside primer binding "
            "regions. It must contain at least Chr and position."
        ),
    )
    args = parser.parse_args()

    sequences = {}
    with open(args.genome_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            chrom_key = normalize_chrom_name(record.id)
            sequences[chrom_key] = str(record.seq).upper()

    marker_info = load_marker_table(args.marker_csv)
    background_positions = build_background_positions(
        marker_info,
        background_csv=args.background_csv,
        background_vcf=args.background_vcf,
    )

    md = MarkerDesign(background_positions=background_positions)
    data = []
    skipped = 0
    skip_reasons = defaultdict(int)
    ref_mismatch_count = 0

    for row in tqdm(
        marker_info.itertuples(index=False),
        total=marker_info.shape[0],
        desc="Processing markers",
    ):
        result, reason, alignment_check = md.design_kasp(
            sequences,
            row.Chr,
            int(row.position),
            row.ref,
            row.alt,
        )
        if not alignment_check:
            ref_mismatch_count += 1
        if result is None:
            skipped += 1
            skip_reasons[reason] += 1
            continue

        skip_reasons[reason] += 1
        data.append(
            {
                "Chr": row.Chr,
                "Position": int(row.position),
                "Ref": row.ref,
                "Alt": row.alt,
                **result,
            }
        )

    df = pd.DataFrame(data)
    df.to_csv(args.output_csv, index=False)
    print(f"Designed {len(data)} primer sets from {marker_info.shape[0]} targets.")
    print(f"Skipped {skipped} targets with no acceptable primer pair.")
    print(f"Reference mismatches against FASTA: {ref_mismatch_count}")
    if skip_reasons:
        print("Outcome summary:")
        for reason, count in sorted(skip_reasons.items(), key=lambda item: (-item[1], item[0])):
            print(f"  {reason}: {count}")
    print(f"Output saved to {args.output_csv}")
