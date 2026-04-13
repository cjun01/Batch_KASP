import argparse
import gzip
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from bisect import bisect_right
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from Bio import SeqIO
import pandas as pd

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover
    def tqdm(iterable, **_kwargs):
        return iterable

try:
    import primer3  # type: ignore
except ImportError:  # pragma: no cover
    primer3 = None


FAM = "GAAGGTGACCAAGTTCATGCT"
VIC = "GAAGGTCGGAGTCAACGGATT"
VALID_BASES = {"A", "C", "G", "T"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")
CONCISE_OUTPUT_COLUMNS: Sequence[str] = (
    "Chr",
    "Position",
    "Ref",
    "Alt",
    "FASTA_base",
    "FASTA_allele_annotation",
    "TASSEL_major_allele_as_REF",
    "Orientation",
    "A1_Primer",
    "A2_Primer",
    "Common_Primer",
    "A1_Core",
    "A2_Core",
    "Common_Core",
    "Allele_specific_span_start",
    "Allele_specific_span_end",
    "Common_primer_span_start",
    "Common_primer_span_end",
    "A1_core_tm",
    "A2_core_tm",
    "Common_core_tm",
    "A1_core_gc",
    "A2_core_gc",
    "Common_core_gc",
    "Product_size",
    "passed",
    "Exclusion_reason",
    "worst_struct_tm",
    "max_tm_delta",
    "Local_score",
    "Ranking_mode",
    "Rank",
    "Design_status",
)
CONCISE_BLAST_OUTPUT_COLUMNS: Sequence[str] = (
    "Total_exact_full_length_off_target_hits",
    "Total_off_target_3prime_hits",
    "Worst_offtarget_3prime_bitscore",
)
TALL_OUTPUT_COLUMNS: Sequence[str] = (
    "Primer_ID",
    "Primer_Sequence",
)
TALL_PRIMER_COLUMNS: Sequence[str] = (
    "A1_Primer",
    "A2_Primer",
    "Common_Primer",
    "Reverse_Primer",
)
FAILED_OUTPUT_COLUMNS: Sequence[str] = (
    "Chr",
    "Position",
    "Ref",
    "Alt",
    "FASTA_base",
    "FASTA_allele_annotation",
    "TASSEL_major_allele_as_REF",
    "Failure_reason",
    "Failure_reason_counts",
    "PASS_candidates_found",
    "FALLBACK_candidates_found",
)


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def normalize_chrom_name(name) -> str:
    value = str(name).strip()
    lowered = value.lower()
    if lowered.startswith("chr"):
        return lowered[3:]
    return lowered


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def is_clean_sequence(seq: str) -> bool:
    return bool(seq) and set(seq).issubset(VALID_BASES)


def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    seq = seq.upper()
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


def strip_known_sequence_suffix(path_str: str) -> str:
    path = Path(path_str)
    name = path.name
    lower = name.lower()
    if lower.endswith(".gz"):
        name = name[:-3]
        lower = lower[:-3]
    for suffix in (".fasta", ".fa", ".fna", ".fas", ".ffn", ".faa"):
        if lower.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return str(path.with_name(name))


def looks_like_fasta_path(path_str: str) -> bool:
    lower = path_str.lower()
    return lower.endswith((".fa", ".fasta", ".fna", ".fas", ".ffn", ".fa.gz", ".fasta.gz", ".fna.gz", ".fas.gz", ".ffn.gz"))


def blast_db_exists(prefix: str) -> bool:
    known_suffixes = [".nhr", ".nin", ".nsq", ".ndb", ".not", ".ntf", ".nto", ".nog", ".nos"]
    return any(Path(prefix + suffix).exists() for suffix in known_suffixes)


def ensure_blast_db(
    blast_reference: str,
    makeblastdb_path: str = "makeblastdb",
    blast_db_dir: Optional[str] = None,
) -> Tuple[str, bool]:
    if blast_db_exists(blast_reference):
        return blast_reference, False

    reference_path = Path(blast_reference)
    if not reference_path.exists() and not looks_like_fasta_path(blast_reference):
        raise FileNotFoundError(
            f"Could not find BLAST database prefix '{blast_reference}' and it does not look like a FASTA file path."
        )

    if shutil.which(makeblastdb_path) is None:
        raise FileNotFoundError(
            f"Could not find makeblastdb executable '{makeblastdb_path}' in PATH."
        )

    if not reference_path.exists():
        raise FileNotFoundError(f"FASTA file for BLAST database creation was not found: {blast_reference}")

    if blast_db_dir:
        output_dir = Path(blast_db_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        db_prefix = output_dir / Path(strip_known_sequence_suffix(reference_path.name)).name
    else:
        db_prefix = Path(strip_known_sequence_suffix(str(reference_path)))

    db_prefix_str = str(db_prefix)
    if blast_db_exists(db_prefix_str):
        return db_prefix_str, False

    cmd = [
        makeblastdb_path,
        "-in",
        str(reference_path),
        "-dbtype",
        "nucl",
        "-parse_seqids",
        "-out",
        db_prefix_str,
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    return db_prefix_str, True


class BackgroundIntervals:
    def __init__(self):
        self._intervals = defaultdict(list)
        self._starts = {}

    def add(self, chrom: str, start: int, span: int = 1):
        start = int(start)
        span = max(1, int(span))
        end = start + span - 1
        self._intervals[normalize_chrom_name(chrom)].append((start, end))

    def finalize(self):
        for chrom, intervals in self._intervals.items():
            intervals.sort(key=lambda item: (item[0], item[1]))
            merged = []
            for start, end in intervals:
                if not merged or merged[-1][1] < start - 1:
                    merged.append((start, end))
                else:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            self._intervals[chrom] = merged
            self._starts[chrom] = [interval[0] for interval in merged]

    def has_overlap(
        self,
        chrom: str,
        start: int,
        end: int,
        ignore_positions: Optional[Set[int]] = None,
    ) -> bool:
        chrom_key = normalize_chrom_name(chrom)
        intervals = self._intervals.get(chrom_key)
        if not intervals:
            return False

        starts = self._starts[chrom_key]
        idx = bisect_right(starts, end)
        if idx == 0:
            return False

        ignored = set(ignore_positions or set())
        for i in range(idx - 1, -1, -1):
            iv_start, iv_end = intervals[i]
            if iv_end < start:
                break
            overlap_start = max(start, iv_start)
            overlap_end = min(end, iv_end)
            if overlap_start > overlap_end:
                continue
            if not ignored:
                return True
            overlap_len = overlap_end - overlap_start + 1
            ignored_len = sum(1 for pos in ignored if overlap_start <= pos <= overlap_end)
            if ignored_len < overlap_len:
                return True
        return False


def normalize_marker_columns(df: pd.DataFrame) -> pd.DataFrame:
    column_map = {column.lower(): column for column in df.columns}
    required = {"chr", "position", "ref", "alt"}
    missing = required - set(column_map)
    if missing:
        raise ValueError(
            "Marker CSV must contain the columns Chr, position, ref, and alt. "
            f"Missing: {', '.join(sorted(missing))}"
        )

    normalized = pd.DataFrame(
        {
            "Chr": df[column_map["chr"]].astype(str),
            "Position": pd.to_numeric(df[column_map["position"]], errors="raise").astype(int),
            "Ref": df[column_map["ref"]].astype(str).str.upper(),
            "Alt": df[column_map["alt"]].astype(str).str.upper(),
        }
    )

    valid_mask = (
        normalized["Ref"].str.len().eq(1)
        & normalized["Alt"].str.len().eq(1)
        & normalized["Ref"].isin(VALID_BASES)
        & normalized["Alt"].isin(VALID_BASES)
        & normalized["Ref"].ne(normalized["Alt"])
    )

    skipped = int((~valid_mask).sum())
    if skipped:
        print(f"Skipping {skipped} non-SNP or invalid markers from the marker CSV.")

    normalized = normalized.loc[valid_mask].drop_duplicates().sort_values(["Chr", "Position"])
    return normalized.reset_index(drop=True)


def load_marker_table(path: str) -> pd.DataFrame:
    return normalize_marker_columns(pd.read_csv(path))


def load_background_from_csv(path: str, background: BackgroundIntervals):
    df = pd.read_csv(path)
    column_map = {column.lower(): column for column in df.columns}
    required = {"chr", "position"}
    missing = required - set(column_map)
    if missing:
        raise ValueError(
            "Background CSV must contain at least the columns Chr and position. "
            f"Missing: {', '.join(sorted(missing))}"
        )

    ref_column = column_map.get("ref")
    positions = pd.to_numeric(df[column_map["position"]], errors="raise").astype(int)
    chromosomes = df[column_map["chr"]].astype(str)

    for index, chrom in chromosomes.items():
        start = int(positions.loc[index])
        if ref_column is not None:
            ref_value = str(df.at[index, ref_column]).upper()
            span = len(ref_value) if ref_value and ref_value != "NAN" else 1
        else:
            span = 1
        background.add(chrom, start, span)


def load_background_from_vcf(path: str, background: BackgroundIntervals):
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
            background.add(chrom, start, span)


def build_background_intervals(
    targets: pd.DataFrame,
    background_csv: Optional[str] = None,
    background_vcf: Optional[str] = None,
) -> BackgroundIntervals:
    background = BackgroundIntervals()

    for row in targets.itertuples(index=False):
        background.add(row.Chr, int(row.Position), len(row.Ref))

    if background_csv:
        load_background_from_csv(background_csv, background)
    if background_vcf:
        load_background_from_vcf(background_vcf, background)

    background.finalize()
    return background


class Primer3Thermo:
    def __init__(
        self,
        mv_conc: float = 50.0,
        dv_conc: float = 1.5,
        dntp_conc: float = 0.6,
        dna_conc: float = 50.0,
        temp_c: float = 37.0,
        max_loop: int = 30,
    ):
        if primer3 is None:
            raise ImportError(
                "primer3-py is required. Install it with: pip install primer3-py"
            )
        self.mv_conc = mv_conc
        self.dv_conc = dv_conc
        self.dntp_conc = dntp_conc
        self.dna_conc = dna_conc
        self.temp_c = temp_c
        self.max_loop = max_loop

    def tm(self, seq: str) -> float:
        return float(
            primer3.calc_tm(
                seq,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
            )
        )

    def hairpin(self, seq: str):
        return primer3.calc_hairpin(
            seq,
            mv_conc=self.mv_conc,
            dv_conc=self.dv_conc,
            dntp_conc=self.dntp_conc,
            dna_conc=self.dna_conc,
            temp_c=self.temp_c,
            max_loop=self.max_loop,
        )

    def homodimer(self, seq: str):
        return primer3.calc_homodimer(
            seq,
            mv_conc=self.mv_conc,
            dv_conc=self.dv_conc,
            dntp_conc=self.dntp_conc,
            dna_conc=self.dna_conc,
            temp_c=self.temp_c,
            max_loop=self.max_loop,
        )

    def heterodimer(self, seq1: str, seq2: str):
        return primer3.calc_heterodimer(
            seq1,
            seq2,
            mv_conc=self.mv_conc,
            dv_conc=self.dv_conc,
            dntp_conc=self.dntp_conc,
            dna_conc=self.dna_conc,
            temp_c=self.temp_c,
            max_loop=self.max_loop,
        )


class BlastPrimerChecker:
    def __init__(
        self,
        blast_db: str,
        blastn_path: str = "blastn",
        task: str = "blastn-short",
        word_size: int = 7,
        max_target_seqs: int = 200,
        num_threads: int = 1,
    ):
        if shutil.which(blastn_path) is None:
            raise FileNotFoundError(
                f"Could not find BLAST executable '{blastn_path}' in PATH."
            )
        self.blast_db = blast_db
        self.blastn_path = blastn_path
        self.task = task
        self.word_size = word_size
        self.max_target_seqs = max_target_seqs
        self.num_threads = num_threads

    @staticmethod
    def _is_on_target(hit: dict, expected_chrom: str, expected_start: int, expected_end: int) -> bool:
        if normalize_chrom_name(hit["sseqid"]) != normalize_chrom_name(expected_chrom):
            return False
        hit_start = min(hit["sstart"], hit["send"])
        hit_end = max(hit["sstart"], hit["send"])
        overlap = max(0, min(hit_end, expected_end) - max(hit_start, expected_start) + 1)
        return overlap >= max(8, hit["effective_len"] - 2)

    def screen_primers(self, primer_queries: List[dict]) -> Dict[str, dict]:
        if not primer_queries:
            return {}

        with tempfile.TemporaryDirectory(prefix="kasp_blast_") as tmpdir:
            query_fa = f"{tmpdir}/queries.fa"
            out_tsv = f"{tmpdir}/hits.tsv"
            with open(query_fa, "w") as handle:
                for query in primer_queries:
                    handle.write(f">{query['qid']}\n{query['seq']}\n")

            cmd = [
                self.blastn_path,
                "-task", self.task,
                "-db", self.blast_db,
                "-query", query_fa,
                "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen",
                "-dust", "no",
                "-soft_masking", "false",
                "-word_size", str(self.word_size),
                "-max_target_seqs", str(self.max_target_seqs),
                "-num_threads", str(self.num_threads),
                "-out", out_tsv,
            ]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)

            by_qid = defaultdict(list)
            with open(out_tsv, "r") as handle:
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) != 13:
                        continue
                    qid = fields[0]
                    qlen = int(fields[12])
                    qstart = int(fields[6])
                    qend = int(fields[7])
                    align_len = int(fields[3])
                    effective_len = qend - qstart + 1
                    hit = {
                        "qseqid": qid,
                        "sseqid": fields[1],
                        "pident": float(fields[2]),
                        "length": align_len,
                        "mismatch": int(fields[4]),
                        "gapopen": int(fields[5]),
                        "qstart": qstart,
                        "qend": qend,
                        "sstart": int(fields[8]),
                        "send": int(fields[9]),
                        "evalue": float(fields[10]),
                        "bitscore": float(fields[11]),
                        "qlen": qlen,
                        "effective_len": effective_len,
                        "covers_3prime_end": qend == qlen,
                        "near_full_length": effective_len >= max(qlen - 2, 10),
                        "exact_full_length": qstart == 1 and qend == qlen and int(fields[4]) == 0 and int(fields[5]) == 0,
                    }
                    by_qid[qid].append(hit)

        summaries = {}
        for query in primer_queries:
            qid = query["qid"]
            hits = by_qid.get(qid, [])
            on_target = []
            off_target = []
            for hit in hits:
                if self._is_on_target(hit, query["chrom"], query["start"], query["end"]):
                    on_target.append(hit)
                else:
                    off_target.append(hit)

            off_target_3prime = [
                hit for hit in off_target
                if hit["covers_3prime_end"] and hit["near_full_length"]
            ]
            exact_full_off_target = [
                hit for hit in off_target
                if hit["exact_full_length"]
            ]
            best_offtarget_bitscore = max((hit["bitscore"] for hit in off_target_3prime), default=0.0)
            best_offtarget_identity = max((hit["pident"] for hit in off_target_3prime), default=0.0)
            summaries[qid] = {
                "total_hits": len(hits),
                "on_target_hits": len(on_target),
                "off_target_hits": len(off_target),
                "off_target_3prime_hits": len(off_target_3prime),
                "exact_full_length_off_target_hits": len(exact_full_off_target),
                "best_offtarget_3prime_bitscore": best_offtarget_bitscore,
                "best_offtarget_3prime_identity": best_offtarget_identity,
            }
        return summaries


class KASPDesigner:
    def __init__(
        self,
        background: BackgroundIntervals,
        thermo: Primer3Thermo,
        allele_primer_min: int = 18,
        allele_primer_max: int = 30,
        common_primer_min: int = 18,
        common_primer_max: int = 30,
        tm_min: float = 57.0,
        tm_max: float = 65.0,
        tm_tolerance: float = 2.0,
        product_min: int = 60,
        product_max: int = 150,
        scan_window: int = 200,
        min_gc: float = 30.0,
        max_gc: float = 70.0,
        max_thermo_tm: float = 45.0,
        max_candidates_per_orientation: int = 200,
    ):
        self.background = background
        self.thermo = thermo
        self.allele_primer_min = allele_primer_min
        self.allele_primer_max = allele_primer_max
        self.common_primer_min = common_primer_min
        self.common_primer_max = common_primer_max
        self.tm_min = tm_min
        self.tm_max = tm_max
        self.tm_tolerance = tm_tolerance
        self.product_min = product_min
        self.product_max = product_max
        self.scan_window = scan_window
        self.min_gc = min_gc
        self.max_gc = max_gc
        self.max_thermo_tm = max_thermo_tm
        self.max_candidates_per_orientation = max_candidates_per_orientation
        self.ideal_product_size = (self.product_min + self.product_max) / 2.0
        self.ideal_tm = (self.tm_min + self.tm_max) / 2.0
        self.ideal_gc = 50.0

    def classify_marker(self, sequence: str, pos0: int, ref: str, alt: str):
        fasta_base = sequence[pos0]
        ref_matches = fasta_base == ref
        alt_matches = fasta_base == alt
        if ref_matches:
            return True, fasta_base, "REF matches FASTA", "No"
        if alt_matches:
            return True, fasta_base, "ALT matches FASTA; FASTA carries ALT allele", "Yes"
        return False, fasta_base, "Neither REF nor ALT matches FASTA", "Unknown"

    def candidate_exclusion_reason(self, product_size: int, worst_struct_tm: float) -> str:
        reasons = []

        if product_size < self.product_min:
            reasons.append("product_size_below_min")
        elif product_size > self.product_max:
            reasons.append("product_size_above_max")

        if worst_struct_tm > self.max_thermo_tm:
            reasons.append("worst_struct_tm_above_max")

        return ";".join(reasons)
    def common_reverse_candidates(
        self,
        chrom: str,
        right_flank: str,
        right_flank_start_1b: int,
        target_tm: float,
        target_pos_1b: int,
    ) -> Iterable[dict]:
        yielded = 0
        for start_offset in range(0, len(right_flank) - self.common_primer_min + 1):
            for primer_len in range(self.common_primer_min, self.common_primer_max + 1):
                end_offset = start_offset + primer_len
                if end_offset > len(right_flank):
                    break
                binding_sequence = right_flank[start_offset:end_offset]
                if not is_clean_sequence(binding_sequence):
                    continue
                core = reverse_complement(binding_sequence)
                tm = self.thermo.tm(core)
                gc = gc_content(core)
                if not (target_tm - self.tm_tolerance <= tm <= target_tm + self.tm_tolerance):
                    continue
                if not (self.min_gc <= gc <= self.max_gc):
                    continue
                start_1b = right_flank_start_1b + start_offset
                end_1b = right_flank_start_1b + end_offset - 1
                if self.background.has_overlap(
                    chrom,
                    start_1b,
                    end_1b,
                    ignore_positions={target_pos_1b},
                ):
                    continue
                yielded += 1
                yield {
                    "core": core,
                    "tm": tm,
                    "gc": gc,
                    "start": start_1b,
                    "end": end_1b,
                    "direction": "reverse",
                }
                if yielded >= self.max_candidates_per_orientation:
                    return

    def common_forward_candidates(
        self,
        chrom: str,
        left_flank: str,
        left_flank_start_1b: int,
        target_tm: float,
        target_pos_1b: int,
    ) -> Iterable[dict]:
        yielded = 0
        for start_offset in range(0, len(left_flank) - self.common_primer_min + 1):
            for primer_len in range(self.common_primer_min, self.common_primer_max + 1):
                end_offset = start_offset + primer_len
                if end_offset > len(left_flank):
                    break
                core = left_flank[start_offset:end_offset]
                if not is_clean_sequence(core):
                    continue
                tm = self.thermo.tm(core)
                gc = gc_content(core)
                if not (target_tm - self.tm_tolerance <= tm <= target_tm + self.tm_tolerance):
                    continue
                if not (self.min_gc <= gc <= self.max_gc):
                    continue
                start_1b = left_flank_start_1b + start_offset
                end_1b = left_flank_start_1b + end_offset - 1
                if self.background.has_overlap(
                    chrom,
                    start_1b,
                    end_1b,
                    ignore_positions={target_pos_1b},
                ):
                    continue
                yielded += 1
                yield {
                    "core": core,
                    "tm": tm,
                    "gc": gc,
                    "start": start_1b,
                    "end": end_1b,
                    "direction": "forward",
                }
                if yielded >= self.max_candidates_per_orientation:
                    return

    def evaluate_pair(
        self,
        a1_full: str,
        a2_full: str,
        common_full: str,
        a1_core_tm: float,
        a2_core_tm: float,
        common_core_tm: float,
        product_size: int,
    ) -> dict:
        a1_hairpin = self.thermo.hairpin(a1_full)
        a2_hairpin = self.thermo.hairpin(a2_full)
        common_hairpin = self.thermo.hairpin(common_full)

        a1_homodimer = self.thermo.homodimer(a1_full)
        a2_homodimer = self.thermo.homodimer(a2_full)
        common_homodimer = self.thermo.homodimer(common_full)

        a1_common_heterodimer = self.thermo.heterodimer(a1_full, common_full)
        a2_common_heterodimer = self.thermo.heterodimer(a2_full, common_full)
        a1_a2_heterodimer = self.thermo.heterodimer(a1_full, a2_full)

        structural_tms = [
            float(a1_hairpin.tm),
            float(a2_hairpin.tm),
            float(common_hairpin.tm),
            float(a1_homodimer.tm),
            float(a2_homodimer.tm),
            float(common_homodimer.tm),
            float(a1_common_heterodimer.tm),
            float(a2_common_heterodimer.tm),
            float(a1_a2_heterodimer.tm),
        ]
        worst_struct_tm = max(structural_tms)
        max_tm_delta = max(
            abs(common_core_tm - a1_core_tm),
            abs(common_core_tm - a2_core_tm),
            abs(a1_core_tm - a2_core_tm),
        )

        exclusion_reason = self.candidate_exclusion_reason(product_size, worst_struct_tm)
        passed = exclusion_reason == ""

        return {
            "passed": passed,
            "Exclusion_reason": exclusion_reason,
            "worst_struct_tm": worst_struct_tm,
            "max_tm_delta": max_tm_delta,
            "A1_hairpin_tm": float(a1_hairpin.tm),
            "A2_hairpin_tm": float(a2_hairpin.tm),
            "Common_hairpin_tm": float(common_hairpin.tm),
            "A1_homodimer_tm": float(a1_homodimer.tm),
            "A2_homodimer_tm": float(a2_homodimer.tm),
            "Common_homodimer_tm": float(common_homodimer.tm),
            "A1_Common_heterodimer_tm": float(a1_common_heterodimer.tm),
            "A2_Common_heterodimer_tm": float(a2_common_heterodimer.tm),
            "A1_A2_heterodimer_tm": float(a1_a2_heterodimer.tm),
        }
    def local_score_tuple(self, result: dict) -> Tuple[float, float, float, float, float]:
        gc_penalty = (
            abs(result["A1_core_gc"] - self.ideal_gc)
            + abs(result["A2_core_gc"] - self.ideal_gc)
            + abs(result["Common_core_gc"] - self.ideal_gc)
        )
        tm_center_penalty = (
            abs(result["A1_core_tm"] - self.ideal_tm)
            + abs(result["A2_core_tm"] - self.ideal_tm)
            + abs(result["Common_core_tm"] - self.ideal_tm)
        )
        return (
            result["worst_struct_tm"],
            result["max_tm_delta"],
            abs(result["Product_size"] - self.ideal_product_size),
            gc_penalty,
            tm_center_penalty,
        )

    def local_score_value(self, result: dict) -> float:
        gc_penalty = (
            abs(result["A1_core_gc"] - self.ideal_gc)
            + abs(result["A2_core_gc"] - self.ideal_gc)
            + abs(result["Common_core_gc"] - self.ideal_gc)
        )
        tm_center_penalty = (
            abs(result["A1_core_tm"] - self.ideal_tm)
            + abs(result["A2_core_tm"] - self.ideal_tm)
            + abs(result["Common_core_tm"] - self.ideal_tm)
        )
        return (
            result["worst_struct_tm"] * 10.0
            + result["max_tm_delta"] * 5.0
            + abs(result["Product_size"] - self.ideal_product_size) * 0.5
            + gc_penalty * 0.25
            + tm_center_penalty * 0.5
        )

    def _finalize_candidate(self, result: dict):
        result["Local_score"] = round(self.local_score_value(result), 4)
        result["Local_score_tuple"] = "|".join(f"{x:.4f}" for x in self.local_score_tuple(result))
        return result

    def design_marker_candidates(self, sequences: Dict[str, str], chrom: str, position: int, ref: str, alt: str):
        sequence = sequences.get(normalize_chrom_name(chrom))
        if sequence is None:
            return {"status": "missing_chromosome_in_fasta", "pass_candidates": [], "fallback_candidates": [], "info": None}

        pos1 = int(position)
        pos0 = pos1 - 1
        if pos0 < 0 or pos0 >= len(sequence):
            return {"status": "position_out_of_range", "pass_candidates": [], "fallback_candidates": [], "info": None}

        ref = ref.upper()
        alt = alt.upper()
        alignment_ok, fasta_base, allele_note, tassel_major_as_ref = self.classify_marker(sequence, pos0, ref, alt)
        base_info = {
            "Chr": chrom,
            "Position": pos1,
            "Ref": ref,
            "Alt": alt,
            "FASTA_base": fasta_base,
            "FASTA_allele_annotation": allele_note,
            "TASSEL_major_allele_as_REF": tassel_major_as_ref,
        }
        if not alignment_ok:
            return {
                "status": "ref_alt_do_not_match_fasta",
                "pass_candidates": [],
                "fallback_candidates": [],
                "info": base_info,
            }

        pass_candidates = []
        fallback_candidates = []
        failure_reasons = defaultdict(int)
        seen_keys = set()

        right_flank = sequence[pos0 + 1: pos0 + 1 + self.scan_window]
        left_flank_start0 = max(0, pos0 - self.scan_window)
        left_flank = sequence[left_flank_start0:pos0]

        # Orientation 1: allele-specific forward, common reverse.
        for allele_len in range(self.allele_primer_min, self.allele_primer_max + 1):
            start0 = pos0 - (allele_len - 1)
            if start0 < 0:
                continue
            template = sequence[start0:pos0 + 1]
            if not is_clean_sequence(template[:-1]):
                failure_reasons["ambiguous_base_in_forward_allele_specific_region"] += 1
                continue
            ref_core = template[:-1] + ref
            alt_core = template[:-1] + alt
            if not is_clean_sequence(ref_core) or not is_clean_sequence(alt_core):
                failure_reasons["invalid_forward_allele_specific_sequence"] += 1
                continue
            ref_tm = self.thermo.tm(ref_core)
            alt_tm = self.thermo.tm(alt_core)
            gc = gc_content(ref_core)
            if not (self.tm_min <= ref_tm <= self.tm_max and self.tm_min <= alt_tm <= self.tm_max):
                failure_reasons["forward_allele_specific_tm"] += 1
                continue
            if not (self.min_gc <= gc <= self.max_gc):
                failure_reasons["forward_allele_specific_gc"] += 1
                continue
            if self.background.has_overlap(chrom, start0 + 1, pos0 + 1, ignore_positions={pos1}):
                failure_reasons["forward_allele_specific_overlaps_variant"] += 1
                continue

            a1_full = FAM + ref_core
            a2_full = VIC + alt_core
            target_tm = (ref_tm + alt_tm) / 2.0
            for common in self.common_reverse_candidates(chrom, right_flank, pos0 + 2, target_tm, pos1):
                key = ("F", ref_core, alt_core, common["core"])
                if key in seen_keys:
                    continue
                seen_keys.add(key)
                product_size = common["end"] - (start0 + 1) + 1
                pair_eval = self.evaluate_pair(
                    a1_full,
                    a2_full,
                    common["core"],
                    ref_tm,
                    alt_tm,
                    common["tm"],
                    product_size,
                )
                result = {
                    **base_info,
                    "Orientation": "allele_specific_forward_common_reverse",
                    "A1_Primer": a1_full,
                    "A2_Primer": a2_full,
                    "Common_Primer": common["core"],
                    "A1_Core": ref_core,
                    "A2_Core": alt_core,
                    "Common_Core": common["core"],
                    "Allele_specific_span_start": start0 + 1,
                    "Allele_specific_span_end": pos0 + 1,
                    "Common_primer_span_start": common["start"],
                    "Common_primer_span_end": common["end"],
                    "A1_core_tm": ref_tm,
                    "A2_core_tm": alt_tm,
                    "Common_core_tm": common["tm"],
                    "A1_core_gc": gc_content(ref_core),
                    "A2_core_gc": gc_content(alt_core),
                    "Common_core_gc": common["gc"],
                    "Product_size": product_size,
                    **pair_eval,
                }
                result = self._finalize_candidate(result)
                if pair_eval["passed"]:
                    pass_candidates.append(result)
                else:
                    fallback_candidates.append(result)

        # Orientation 2: allele-specific reverse, common forward.
        for allele_len in range(self.allele_primer_min, self.allele_primer_max + 1):
            end0 = pos0 + allele_len
            if end0 > len(sequence):
                continue
            template = sequence[pos0:end0]
            if len(template) != allele_len:
                continue
            if not is_clean_sequence(template[1:]):
                failure_reasons["ambiguous_base_in_reverse_allele_specific_region"] += 1
                continue
            ref_core = reverse_complement(ref + template[1:])
            alt_core = reverse_complement(alt + template[1:])
            if not is_clean_sequence(ref_core) or not is_clean_sequence(alt_core):
                failure_reasons["invalid_reverse_allele_specific_sequence"] += 1
                continue
            ref_tm = self.thermo.tm(ref_core)
            alt_tm = self.thermo.tm(alt_core)
            gc = gc_content(ref_core)
            if not (self.tm_min <= ref_tm <= self.tm_max and self.tm_min <= alt_tm <= self.tm_max):
                failure_reasons["reverse_allele_specific_tm"] += 1
                continue
            if not (self.min_gc <= gc <= self.max_gc):
                failure_reasons["reverse_allele_specific_gc"] += 1
                continue
            if self.background.has_overlap(chrom, pos1, end0, ignore_positions={pos1}):
                failure_reasons["reverse_allele_specific_overlaps_variant"] += 1
                continue

            a1_full = FAM + ref_core
            a2_full = VIC + alt_core
            target_tm = (ref_tm + alt_tm) / 2.0
            for common in self.common_forward_candidates(chrom, left_flank, left_flank_start0 + 1, target_tm, pos1):
                key = ("R", ref_core, alt_core, common["core"])
                if key in seen_keys:
                    continue
                seen_keys.add(key)
                product_size = end0 - common["start"] + 1
                pair_eval = self.evaluate_pair(
                    a1_full,
                    a2_full,
                    common["core"],
                    ref_tm,
                    alt_tm,
                    common["tm"],
                    product_size,
                )
                result = {
                    **base_info,
                    "Orientation": "allele_specific_reverse_common_forward",
                    "A1_Primer": a1_full,
                    "A2_Primer": a2_full,
                    "Common_Primer": common["core"],
                    "A1_Core": ref_core,
                    "A2_Core": alt_core,
                    "Common_Core": common["core"],
                    "Allele_specific_span_start": pos1,
                    "Allele_specific_span_end": end0,
                    "Common_primer_span_start": common["start"],
                    "Common_primer_span_end": common["end"],
                    "A1_core_tm": ref_tm,
                    "A2_core_tm": alt_tm,
                    "Common_core_tm": common["tm"],
                    "A1_core_gc": gc_content(ref_core),
                    "A2_core_gc": gc_content(alt_core),
                    "Common_core_gc": common["gc"],
                    "Product_size": product_size,
                    **pair_eval,
                }
                result = self._finalize_candidate(result)
                if pair_eval["passed"]:
                    pass_candidates.append(result)
                else:
                    fallback_candidates.append(result)

        pass_candidates.sort(key=self.local_score_tuple)
        fallback_candidates.sort(key=self.local_score_tuple)
        if pass_candidates:
            status = "designed"
        elif fallback_candidates:
            status = "thermo_fallback"
        elif failure_reasons:
            status = max(failure_reasons.items(), key=lambda item: item[1])[0]
        else:
            status = "no_candidate_found"

        return {
            "status": status,
            "pass_candidates": pass_candidates,
            "fallback_candidates": fallback_candidates,
            "info": base_info,
            "failure_reasons": dict(failure_reasons),
        }


def load_sequences(genome_path: str) -> Dict[str, str]:
    sequences = {}
    with open_text(genome_path) as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            sequences[normalize_chrom_name(record.id)] = str(record.seq).upper()
    return sequences


def build_blast_queries(candidates: List[dict]) -> List[dict]:
    queries = []
    for idx, candidate in enumerate(candidates):
        query_base = f"cand{idx}"
        queries.append(
            {
                "qid": f"{query_base}|A1",
                "seq": candidate["A1_Core"],
                "chrom": candidate["Chr"],
                "start": int(candidate["Allele_specific_span_start"]),
                "end": int(candidate["Allele_specific_span_end"]),
                "role": "A1",
            }
        )
        queries.append(
            {
                "qid": f"{query_base}|A2",
                "seq": candidate["A2_Core"],
                "chrom": candidate["Chr"],
                "start": int(candidate["Allele_specific_span_start"]),
                "end": int(candidate["Allele_specific_span_end"]),
                "role": "A2",
            }
        )
        queries.append(
            {
                "qid": f"{query_base}|COMMON",
                "seq": candidate["Common_Core"],
                "chrom": candidate["Chr"],
                "start": int(candidate["Common_primer_span_start"]),
                "end": int(candidate["Common_primer_span_end"]),
                "role": "COMMON",
            }
        )
    return queries


def annotate_offtarget_metrics(candidates: List[dict], checker: BlastPrimerChecker):
    if not candidates:
        return
    queries = build_blast_queries(candidates)
    summaries = checker.screen_primers(queries)

    for idx, candidate in enumerate(candidates):
        a1 = summaries.get(f"cand{idx}|A1", {})
        a2 = summaries.get(f"cand{idx}|A2", {})
        common = summaries.get(f"cand{idx}|COMMON", {})
        candidate["A1_off_target_3prime_hits"] = a1.get("off_target_3prime_hits", 0)
        candidate["A2_off_target_3prime_hits"] = a2.get("off_target_3prime_hits", 0)
        candidate["Common_off_target_3prime_hits"] = common.get("off_target_3prime_hits", 0)
        candidate["A1_exact_full_length_off_target_hits"] = a1.get("exact_full_length_off_target_hits", 0)
        candidate["A2_exact_full_length_off_target_hits"] = a2.get("exact_full_length_off_target_hits", 0)
        candidate["Common_exact_full_length_off_target_hits"] = common.get("exact_full_length_off_target_hits", 0)
        candidate["A1_best_offtarget_3prime_bitscore"] = a1.get("best_offtarget_3prime_bitscore", 0.0)
        candidate["A2_best_offtarget_3prime_bitscore"] = a2.get("best_offtarget_3prime_bitscore", 0.0)
        candidate["Common_best_offtarget_3prime_bitscore"] = common.get("best_offtarget_3prime_bitscore", 0.0)
        candidate["Total_off_target_3prime_hits"] = (
            candidate["A1_off_target_3prime_hits"]
            + candidate["A2_off_target_3prime_hits"]
            + candidate["Common_off_target_3prime_hits"]
        )
        candidate["Total_exact_full_length_off_target_hits"] = (
            candidate["A1_exact_full_length_off_target_hits"]
            + candidate["A2_exact_full_length_off_target_hits"]
            + candidate["Common_exact_full_length_off_target_hits"]
        )
        candidate["Worst_offtarget_3prime_bitscore"] = max(
            candidate["A1_best_offtarget_3prime_bitscore"],
            candidate["A2_best_offtarget_3prime_bitscore"],
            candidate["Common_best_offtarget_3prime_bitscore"],
        )


def combined_sort_key(candidate: dict):
    return (
        int(candidate.get("Total_exact_full_length_off_target_hits", 0)),
        int(candidate.get("Total_off_target_3prime_hits", 0)),
        float(candidate.get("Worst_offtarget_3prime_bitscore", 0.0)),
        float(candidate["Local_score"]),
        float(candidate["worst_struct_tm"]),
        float(candidate["max_tm_delta"]),
    )


def common_primer_sites_too_close(candidate_a: dict, candidate_b: dict, min_gap: int) -> bool:
    gap = max(0, int(min_gap))
    start_a = int(candidate_a["Common_primer_span_start"])
    end_a = int(candidate_a["Common_primer_span_end"])
    start_b = int(candidate_b["Common_primer_span_start"])
    end_b = int(candidate_b["Common_primer_span_end"])
    return not (end_a + gap < start_b or end_b + gap < start_a)


def diversify_ranked_candidates(candidates: List[dict], top_n: int, min_common_primer_gap: int) -> List[dict]:
    if top_n <= 0:
        return []

    selected = []
    for candidate in candidates:
        if any(common_primer_sites_too_close(candidate, kept, min_common_primer_gap) for kept in selected):
            continue
        selected.append(candidate)
        if len(selected) >= top_n:
            break
    return selected


def build_output_dataframe(results: List[dict], output_mode: str, blast_enabled: bool) -> pd.DataFrame:
    output_df = pd.DataFrame(results)
    if output_mode == "full":
        return output_df

    selected_columns = list(CONCISE_OUTPUT_COLUMNS)
    if blast_enabled:
        selected_columns.extend(CONCISE_BLAST_OUTPUT_COLUMNS)

    if output_df.empty:
        return pd.DataFrame(columns=selected_columns)

    return output_df.loc[:, [column for column in selected_columns if column in output_df.columns]]


def format_tall_rank(value) -> str:
    if pd.isna(value):
        return "NA"
    if isinstance(value, float) and value.is_integer():
        value = int(value)
    return str(value)


def build_tall_output_dataframe(output_df: pd.DataFrame) -> pd.DataFrame:
    if output_df.empty:
        return pd.DataFrame(columns=TALL_OUTPUT_COLUMNS)

    primer_columns = [column for column in TALL_PRIMER_COLUMNS if column in output_df.columns]
    if not primer_columns:
        return pd.DataFrame(columns=TALL_OUTPUT_COLUMNS)

    tall_rows = []
    for record in output_df.to_dict("records"):
        chrom = str(record.get("Chr", "NA"))
        position = record.get("Position", "NA")
        rank = format_tall_rank(record.get("Rank"))
        for primer_column in primer_columns:
            sequence = record.get(primer_column, "")
            if pd.isna(sequence):
                continue
            sequence = str(sequence).strip()
            if not sequence:
                continue
            tall_rows.append(
                {
                    "Primer_ID": f"{chrom}_{position}_rank{rank}_{primer_column}",
                    "Primer_Sequence": sequence,
                }
            )

    return pd.DataFrame(tall_rows, columns=TALL_OUTPUT_COLUMNS)


def format_failure_reason_counts(failure_reasons: Dict[str, int]) -> str:
    if not failure_reasons:
        return ""
    return ";".join(
        f"{reason}={count}"
        for reason, count in sorted(failure_reasons.items(), key=lambda item: (-item[1], item[0]))
    )


def build_failed_output_path(output_csv: str, failed_output_csv: Optional[str] = None) -> str:
    if failed_output_csv:
        return failed_output_csv

    output_path = Path(output_csv)
    if output_path.suffix:
        return str(output_path.with_name(f"{output_path.stem}_failed{output_path.suffix}"))
    return f"{output_csv}_failed.csv"


def build_failed_output_dataframe(failed_results: List[dict]) -> pd.DataFrame:
    if not failed_results:
        return pd.DataFrame(columns=FAILED_OUTPUT_COLUMNS)

    failed_df = pd.DataFrame(failed_results)
    return failed_df.loc[:, [column for column in FAILED_OUTPUT_COLUMNS if column in failed_df.columns]]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Design KASP primers with strict nearby-variant rejection, Primer3 thermodynamic checks, "
            "multi-candidate ranking, and optional BLAST off-target screening."
        )
    )
    parser.add_argument("-g", "--genome_path", required=True, help="Path to the genome FASTA file.")
    parser.add_argument(
        "-m",
        "--marker_csv",
        required=True,
        help="Path to the marker CSV containing Chr, position, ref, and alt.",
    )
    parser.add_argument("-o", "--output_csv", required=True, help="Path to save the output CSV.")
    parser.add_argument(
        "--background_vcf",
        help="Optional VCF or VCF.gz of background variants to avoid in primer binding regions.",
    )
    parser.add_argument(
        "--background_csv",
        help="Optional background CSV of variants to avoid in primer binding regions.",
    )
    parser.add_argument("--allele_primer_min", type=int, default=18)
    parser.add_argument("--allele_primer_max", type=int, default=30)
    parser.add_argument("--common_primer_min", type=int, default=18)
    parser.add_argument("--common_primer_max", type=int, default=30)
    parser.add_argument("--tm_min", type=float, default=57.0)
    parser.add_argument("--tm_max", type=float, default=65.0)
    parser.add_argument("--tm_tolerance", type=float, default=2.0)
    parser.add_argument("--product_min", type=int, default=60)
    parser.add_argument("--product_max", type=int, default=150)
    parser.add_argument("--scan_window", type=int, default=200)
    parser.add_argument("--min_gc", type=float, default=30.0)
    parser.add_argument("--max_gc", type=float, default=70.0)
    parser.add_argument(
        "--max_thermo_tm",
        type=float,
        default=45.0,
        help="Maximum acceptable Primer3 hairpin/dimer Tm for a PASS design.",
    )
    parser.add_argument(
        "--max_candidates_per_orientation",
        type=int,
        default=200,
        help="Cap on candidate common primers explored per orientation during local search.",
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=5,
        help="Number of ranked PASS assays to keep per SNP in the output.",
    )
    parser.add_argument(
        "--min_common_primer_gap",
        type=int,
        default=5,
        help=(
            "Require reported assays for the same SNP to have common-primer binding sites separated by at "
            "least this many bases. Set to 0 to allow overlapping common primers in the output."
        ),
    )
    parser.add_argument(
        "--include_fallback",
        action="store_true",
        help="If no PASS assays exist for a SNP, include up to top_n FALLBACK assays in the output.",
    )
    parser.add_argument(
        "--preblast_top_n",
        type=int,
        default=10,
        help="When BLAST screening is enabled, only send the best local candidates per SNP to BLAST.",
    )
    parser.add_argument(
        "--output_mode",
        choices=("concise", "full"),
        default="concise",
        help=(
            "Output schema to write. 'concise' keeps the strict selection logic but writes a much smaller CSV; "
            "'full' writes all diagnostic columns."
        ),
    )
    parser.add_argument(
        "--failed_output_csv",
        help=(
            "Optional companion CSV for markers that do not receive any written assay rows. "
            "Defaults to <output_csv stem>_failed.csv when needed."
        ),
    )
    parser.add_argument(
        "--tall_output_csv",
        help=(
            "Optional tall-format companion CSV with one primer per row. "
            "Writes Primer_ID and Primer_Sequence columns, where Primer_ID is "
            "Chr_Position_rank<Rank>_<PrimerColumn>."
        ),
    )
    parser.add_argument(
        "--blast_db",
        help=(
            "Optional BLAST reference. This can be either an existing local BLAST database prefix "
            "or a genome FASTA path; if a FASTA path is given, makeblastdb will be run automatically."
        ),
    )
    parser.add_argument(
        "--blast_use_genome",
        action="store_true",
        help="Use --genome_path as the BLAST reference and automatically build a database if needed.",
    )
    parser.add_argument("--blast_db_dir", help="Optional directory where an auto-built BLAST database should be written.")
    parser.add_argument("--blastn_path", default="blastn", help="Path to blastn executable.")
    parser.add_argument("--makeblastdb_path", default="makeblastdb", help="Path to makeblastdb executable.")
    parser.add_argument("--blast_task", default="blastn-short", help="BLAST task to use for off-target screening.")
    parser.add_argument("--blast_word_size", type=int, default=7)
    parser.add_argument("--blast_max_target_seqs", type=int, default=200)
    parser.add_argument("--blast_num_threads", type=int, default=1)
    parser.add_argument("--mv_conc", type=float, default=50.0, help="Monovalent cation concentration in mM.")
    parser.add_argument("--dv_conc", type=float, default=1.5, help="Divalent cation concentration in mM.")
    parser.add_argument("--dntp_conc", type=float, default=0.6, help="dNTP concentration in mM.")
    parser.add_argument("--dna_conc", type=float, default=50.0, help="DNA concentration in nM.")
    parser.add_argument("--temp_c", type=float, default=37.0, help="Temperature used for Primer3 structure calculations.")
    parser.add_argument("--max_loop", type=int, default=30, help="Maximum loop size for Primer3 structure calculations.")
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    sequences = load_sequences(args.genome_path)
    marker_info = load_marker_table(args.marker_csv)
    background = build_background_intervals(
        marker_info,
        background_csv=args.background_csv,
        background_vcf=args.background_vcf,
    )

    thermo = Primer3Thermo(
        mv_conc=args.mv_conc,
        dv_conc=args.dv_conc,
        dntp_conc=args.dntp_conc,
        dna_conc=args.dna_conc,
        temp_c=args.temp_c,
        max_loop=args.max_loop,
    )
    designer = KASPDesigner(
        background=background,
        thermo=thermo,
        allele_primer_min=args.allele_primer_min,
        allele_primer_max=args.allele_primer_max,
        common_primer_min=args.common_primer_min,
        common_primer_max=args.common_primer_max,
        tm_min=args.tm_min,
        tm_max=args.tm_max,
        tm_tolerance=args.tm_tolerance,
        product_min=args.product_min,
        product_max=args.product_max,
        scan_window=args.scan_window,
        min_gc=args.min_gc,
        max_gc=args.max_gc,
        max_thermo_tm=args.max_thermo_tm,
        max_candidates_per_orientation=args.max_candidates_per_orientation,
    )

    blast_checker = None
    resolved_blast_db = None
    built_blast_db = False
    blast_reference = None
    if args.blast_use_genome:
        blast_reference = args.genome_path
    elif args.blast_db:
        blast_reference = args.blast_db

    if blast_reference:
        resolved_blast_db, built_blast_db = ensure_blast_db(
            blast_reference,
            makeblastdb_path=args.makeblastdb_path,
            blast_db_dir=args.blast_db_dir,
        )
        blast_checker = BlastPrimerChecker(
            blast_db=resolved_blast_db,
            blastn_path=args.blastn_path,
            task=args.blast_task,
            word_size=args.blast_word_size,
            max_target_seqs=args.blast_max_target_seqs,
            num_threads=args.blast_num_threads,
        )

    results = []
    failed_results = []
    outcome_counts = defaultdict(int)
    ref_mismatch_count = 0
    alt_matches_fasta_count = 0

    for row in tqdm(marker_info.itertuples(index=False), total=marker_info.shape[0], desc="Processing markers"):
        marker_result = designer.design_marker_candidates(
            sequences,
            row.Chr,
            int(row.Position),
            row.Ref,
            row.Alt,
        )
        outcome_counts[marker_result["status"]] += 1

        info = marker_result.get("info") or {
            "Chr": row.Chr,
            "Position": int(row.Position),
            "Ref": row.Ref,
            "Alt": row.Alt,
            "FASTA_base": "NA",
            "FASTA_allele_annotation": "NA",
            "TASSEL_major_allele_as_REF": "Unknown",
        }
        if info.get("FASTA_allele_annotation") != "REF matches FASTA":
            ref_mismatch_count += 1
        if info.get("TASSEL_major_allele_as_REF") == "Yes":
            alt_matches_fasta_count += 1

        selected = []
        if marker_result["pass_candidates"]:
            local_top = (
                diversify_ranked_candidates(
                    marker_result["pass_candidates"],
                    max(1, args.preblast_top_n),
                    args.min_common_primer_gap,
                )
                if blast_checker
                else marker_result["pass_candidates"]
            )
            if blast_checker and local_top:
                annotate_offtarget_metrics(local_top, blast_checker)
                local_top.sort(key=combined_sort_key)
                for candidate in local_top:
                    candidate["Ranking_mode"] = "local_plus_blast"
                selected = diversify_ranked_candidates(
                    local_top,
                    max(1, args.top_n),
                    args.min_common_primer_gap,
                )
            else:
                for candidate in marker_result["pass_candidates"]:
                    candidate["Ranking_mode"] = "local_only"
                selected = diversify_ranked_candidates(
                    marker_result["pass_candidates"],
                    max(1, args.top_n),
                    args.min_common_primer_gap,
                )
        elif args.include_fallback and marker_result["fallback_candidates"]:
            fallback_top = (
                diversify_ranked_candidates(
                    marker_result["fallback_candidates"],
                    max(1, args.preblast_top_n),
                    args.min_common_primer_gap,
                )
                if blast_checker
                else marker_result["fallback_candidates"]
            )
            if blast_checker and fallback_top:
                annotate_offtarget_metrics(fallback_top, blast_checker)
                fallback_top.sort(key=combined_sort_key)
                for candidate in fallback_top:
                    candidate["Ranking_mode"] = "fallback_local_plus_blast"
                selected = diversify_ranked_candidates(
                    fallback_top,
                    max(1, args.top_n),
                    args.min_common_primer_gap,
                )
            else:
                for candidate in marker_result["fallback_candidates"]:
                    candidate["Ranking_mode"] = "fallback_local_only"
                selected = diversify_ranked_candidates(
                    marker_result["fallback_candidates"],
                    max(1, args.top_n),
                    args.min_common_primer_gap,
                )
            for candidate in selected:
                candidate["Design_status"] = "FALLBACK"

        if selected:
            selected.sort(key=combined_sort_key if blast_checker else lambda x: x["Local_score"])
            for rank, candidate in enumerate(selected, start=1):
                candidate["Rank"] = rank
                if blast_checker:
                    candidate["Combined_rank_key"] = "|".join(
                        str(x) for x in combined_sort_key(candidate)
                    )
                else:
                    candidate["Combined_rank_key"] = str(candidate["Local_score"])
                if "Design_status" not in candidate:
                    candidate["Design_status"] = "PASS"
                results.append(candidate)
        else:
            failed_results.append(
                {
                    "Chr": info.get("Chr", row.Chr),
                    "Position": int(info.get("Position", row.Position)),
                    "Ref": info.get("Ref", row.Ref),
                    "Alt": info.get("Alt", row.Alt),
                    "FASTA_base": info.get("FASTA_base", "NA"),
                    "FASTA_allele_annotation": info.get("FASTA_allele_annotation", "NA"),
                    "TASSEL_major_allele_as_REF": info.get("TASSEL_major_allele_as_REF", "Unknown"),
                    "Failure_reason": marker_result["status"],
                    "Failure_reason_counts": format_failure_reason_counts(marker_result.get("failure_reasons") or {}),
                    "PASS_candidates_found": len(marker_result["pass_candidates"]),
                    "FALLBACK_candidates_found": len(marker_result["fallback_candidates"]),
                }
            )

    output_df = build_output_dataframe(
        results,
        output_mode=args.output_mode,
        blast_enabled=blast_checker is not None,
    )
    output_df.to_csv(args.output_csv, index=False)

    tall_output_path = None
    if args.tall_output_csv:
        tall_output_path = args.tall_output_csv
        tall_output_df = build_tall_output_dataframe(output_df)
        tall_output_df.to_csv(tall_output_path, index=False)

    failed_output_path = None
    if failed_results or args.failed_output_csv:
        failed_output_path = build_failed_output_path(args.output_csv, args.failed_output_csv)
        failed_output_df = build_failed_output_dataframe(failed_results)
        failed_output_df.to_csv(failed_output_path, index=False)

    print(f"Output rows written: {len(results)}")
    print(f"Output mode: {args.output_mode}")
    print(f"Markers without written output rows: {len(failed_results)}")
    print(f"Reference mismatches against FASTA: {ref_mismatch_count}")
    print(f"ALT matches FASTA and likely treated as REF by TASSEL: {alt_matches_fasta_count}")
    print("Outcome summary:")
    for reason, count in sorted(outcome_counts.items(), key=lambda item: (-item[1], item[0])):
        print(f"  {reason}: {count}")
    if resolved_blast_db:
        print(f"Off-target BLAST screening enabled against database: {resolved_blast_db}")
        if built_blast_db:
            print("BLAST database was auto-built for this run.")
    print(f"Output saved to {args.output_csv}")
    if tall_output_path:
        print(f"Tall-format output saved to {tall_output_path}")
    if failed_output_path:
        print(f"Failed-marker report saved to {failed_output_path}")


if __name__ == "__main__":
    main()
