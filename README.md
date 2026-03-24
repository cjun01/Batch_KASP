# Batch_KASP

`Batch_KASP` is a batch KASP assay design workflow for **biallelic SNP markers**. It takes a **reference genome FASTA** plus a **marker CSV** and returns one or more ranked KASP assay designs per target SNP.

The main design engine is `KASP_design.py`. It:

- reads a reference FASTA or FASTA.gz
- validates and normalizes marker input
- supports both common KASP assay orientations
- rejects primers that overlap nearby known variants
- uses **Primer3 thermodynamic checks** for hairpins and dimer interactions
- keeps multiple valid candidates per marker
- ranks candidates locally
- can optionally re-rank top candidates using **BLAST off-target screening**

The repository also includes:

- `vcf_to_kasp_csv.py` — convert a VCF/VCF.gz into the marker CSV expected by the designer
- `check_markers_against_fasta.py` — verify marker coordinates and alleles against the reference FASTA

---

## Repository contents

```text
Batch_KASP/
├── KASP_design.py
├── vcf_to_kasp_csv.py
├── check_markers_against_fasta.py
├── requirements.txt
└── README.md
```

---

## What the pipeline does

For each SNP marker, `KASP_design.py` performs the following steps:

1. **Loads the reference genome**
   - accepts FASTA or FASTA.gz
   - normalizes chromosome names so `chr1` and `1` are treated as the same chromosome

2. **Reads and cleans the marker table**
   - requires `Chr`, `position`, `ref`, and `alt`
   - treats column names case-insensitively
   - keeps only simple biallelic A/C/G/T SNPs
   - removes exact duplicate marker rows after normalization

3. **Checks allele consistency against the FASTA**
   - `REF matches FASTA`
   - `ALT matches FASTA; FASTA carries ALT allele`
   - `Neither REF nor ALT matches FASTA`

   Markers in the third category are rejected.

4. **Builds a strict variant mask for primer placement**
   - all target SNPs are masked automatically
   - optional background variants from a VCF or CSV are also masked
   - the current target SNP itself is still allowed inside its own allele-specific primer region
   - any other overlapping known variant causes that primer candidate to be rejected

5. **Searches both supported KASP assay layouts**
   - **Orientation 1**: allele-specific forward primers + common reverse primer
   - **Orientation 2**: allele-specific reverse primers + common forward primer

6. **Adds the standard reporter tails**
   - `A1_Primer = FAM tail + REF-specific core`
   - `A2_Primer = VIC tail + ALT-specific core`

   Default tails:

   - FAM tail: `GAAGGTGACCAAGTTCATGCT`
   - VIC tail: `GAAGGTCGGAGTCAACGGATT`

7. **Evaluates thermodynamic quality with Primer3**
   - hairpins
   - homodimers
   - heterodimers
   - Tm balance among primer cores

8. **Ranks local candidates**
   - candidates are scored and sorted using a strict local ranking rule
   - multiple ranked candidates can be retained per SNP

9. **Optionally performs BLAST off-target screening**
   - BLAST is run on `A1_Core`, `A2_Core`, and `Common_Core`
   - top local candidates are re-ranked using off-target hit summaries

---

## Installation

Install the Python requirements:

```bash
pip install -r requirements.txt
pip install primer3-py
pip install tqdm
```

### Notes

- `biopython` and `pandas` come from `requirements.txt`
- `primer3-py` is required for actual thermodynamic primer evaluation
- `tqdm` is optional; without it, the script still runs but without a progress bar
- BLAST screening is optional, but if you use it you also need **NCBI BLAST+** installed:
  - `blastn`
  - `makeblastdb`

Check BLAST availability with:

```bash
blastn -version
makeblastdb -version
```

---

## Input files

## 1) Reference FASTA

Argument:

```bash
-g / --genome_path
```

Accepted formats:

- `genome.fa`
- `genome.fasta`
- `genome.fa.gz`
- `genome.fasta.gz`

Behavior:

- chromosome matching is case-insensitive
- a leading `chr` prefix is ignored during comparison

Examples:

- `chr1` matches `1`
- `Chr05` matches `05`

---

## 2) Marker CSV

Argument:

```bash
-m / --marker_csv
```

Required columns:

- `Chr`
- `position`
- `ref`
- `alt`

Rules:

- column matching is case-insensitive
- `position` must be **1-based**
- `ref` and `alt` must each be exactly one of `A`, `C`, `G`, or `T`
- non-SNP rows, invalid alleles, and duplicate rows are removed

Example:

```csv
Chr,position,ref,alt
1A,123456,A,G
1B,456789,C,T
```

---

## 3) Optional background variants

You can provide known nearby variants that should be avoided in primer-binding regions.

### Background VCF

```bash
--background_vcf input.vcf.gz
```

Accepted:

- `.vcf`
- `.vcf.gz`

The masked span uses the length of the VCF `REF` allele.

### Background CSV

```bash
--background_csv background.csv
```

Must contain at least:

- `Chr`
- `position`

If a `ref` column is present, the masked span uses its length. Otherwise, a single base is masked.

### Important masking rule

The script is intentionally strict:

- the current target SNP is allowed for its own assay
- any **other** masked variant overlapping a primer binding region causes rejection

This means the pipeline prefers **cleaner assays** over maximizing the number of returned designs.

---

## Workflow overview

A recommended workflow is:

1. Convert a VCF to the marker CSV format
2. Check that marker coordinates and REF alleles match the FASTA
3. Run local KASP design
4. Optionally add BLAST off-target screening

---

## Step 1: Convert VCF to marker CSV

Use `vcf_to_kasp_csv.py` to extract SNP markers into the required input format.

Basic conversion:

```bash
python vcf_to_kasp_csv.py \
  -i input.vcf.gz \
  -o target_markers.csv
```

### Example: sample-specific conversion

Keep only markers for a specific sample:

```bash
python vcf_to_kasp_csv.py \
  -i input.vcf.gz \
  -o target_markers.csv \
  -s SampleName
```

Keep only non-reference genotypes for that sample:

```bash
python vcf_to_kasp_csv.py \
  -i input.vcf.gz \
  -o target_markers.csv \
  -s SampleName \
  --genotype nonref
```

Typical output:

```csv
Chr,position,ref,alt
2,10455,G,A
2,20231,T,C
5,99887,C,G
```

Use this output directly as the `-m/--marker_csv` input to `KASP_design.py`.

---

## Step 2: Check markers against the FASTA

Before design, verify that the marker coordinates and `REF` alleles are consistent with the reference.

```bash
python check_markers_against_fasta.py \
  -g genome.fa \
  -m target_markers.csv \
  -o marker_fasta_check.csv
```

This is useful for catching:

- wrong genome version
- coordinate shifts
- chromosome naming mismatches
- cases where `ALT`, not `REF`, matches the FASTA

You should inspect the resulting report before running large design jobs.

---

## Step 3: Run local KASP design

Minimal local design example:

```bash
python KASP_design.py \
  -g genome.fa \
  -m target_markers.csv \
  -o primer_designs.csv
```

By default, the current script now writes a **concise output table**: it keeps the strict internal screening and ranking logic, but only writes a smaller set of decision-useful columns. Use `--output_mode full` if you want the full diagnostic table.

If a marker does not receive any written assay row, the script also writes a companion failed-marker report by default using the output stem plus `_failed.csv`. For example, `primer_designs.csv` is paired with `primer_designs_failed.csv`.

### Example: include a background VCF

This is usually the recommended mode when your markers came from a larger variant set.

```bash
python KASP_design.py \
  -g genome.fa \
  -m target_markers.csv \
  -o primer_designs.csv \
  --background_vcf input.vcf.gz
```

### Example: use a background CSV instead

```bash
python KASP_design.py \
  -g genome.fa \
  -m target_markers.csv \
  -o primer_designs.csv \
  --background_csv background.csv
```

### Example: stricter local filters

```bash
python KASP_design.py \
  -g genome.fa \
  -m target_markers.csv \
  -o primer_designs.csv \
  --background_vcf input.vcf.gz \
  --product_min 60 \
  --product_max 120 \
  --max_thermo_tm 40
```

This configuration narrows the acceptable amplicon size and lowers tolerance for strong unwanted structures.

---

## Step 4: Add BLAST off-target screening

BLAST screening is optional and is mainly useful when:

- the genome is large or repetitive
- duplicated loci are a concern
- you want additional confidence in specificity
- you want BLAST-assisted ranking among already acceptable local candidates

### Option A: use an existing BLAST database

```bash
python KASP_design.py \
  -g genome.fa \
  -m target_markers.csv \
  -o primer_designs.csv \
  --background_vcf input.vcf.gz \
  --blast_db /path/to/blast_db_prefix
```

### Option B: use the genome FASTA directly and auto-build the BLAST database

```bash
python KASP_design.py \
  -g genome.fa \
  -m target_markers.csv \
  -o primer_designs.csv \
  --background_vcf input.vcf.gz \
  --blast_use_genome
```

### Option C: control where the auto-built BLAST database is written

```bash
python KASP_design.py \
  -g genome.fa \
  -m target_markers.csv \
  -o primer_designs.csv \
  --background_vcf input.vcf.gz \
  --blast_use_genome \
  --blast_db_dir blastdb
```

### Option D: use a FASTA path through `--blast_db`

If the code supports FASTA autodetection for `--blast_db`, you can also do:

```bash
python KASP_design.py \
  -g genome.fa \
  -m target_markers.csv \
  -o primer_designs.csv \
  --background_vcf input.vcf.gz \
  --blast_db genome.fa
```

---

## Candidate generation and ranking

The script does **not** stop at the first acceptable assay. Instead, it can keep multiple acceptable candidates per SNP and rank them.

## Local candidate generation

For each SNP:

- both supported orientations are explored
- allele-specific and common primers must satisfy current length, Tm, and GC constraints
- any primer overlapping a masked background variant is rejected, except for the marker's own SNP position
- Primer3 structure checks are computed for each surviving candidate pair
- each candidate gets a local ranking record

### Local ranking logic

Candidates are locally ordered by `Local_score_tuple`, which stores:

1. `worst_struct_tm`
2. `max_tm_delta`
3. distance from the ideal product-size midpoint
4. GC penalty relative to 50%
5. Tm penalty relative to the midpoint of `tm_min` and `tm_max`

The script also writes a scalar `Local_score`, but the primary local ordering is driven by `Local_score_tuple`.

### Why this matters

This means the pipeline can return several acceptable designs for the same SNP and rank them from better to worse, instead of only returning the first design encountered during scanning.

---

## `PASS` vs `FALLBACK`

Each candidate is first assessed locally.

The raw boolean field `passed` is `True` only if both conditions are satisfied:

- `product_min <= Product_size <= product_max`
- `worst_struct_tm <= max_thermo_tm`

Selected output rows use `Design_status`:

- `PASS` — selected from candidates that satisfy the local pass criteria
- `FALLBACK` — selected from weaker candidates only when no `PASS` candidates exist

By default:

- only `PASS` rows are written
- `FALLBACK` rows are omitted

To include fallback rows when a marker has no passing candidate:

```bash
python KASP_design.py \
  -g genome.fa \
  -m target_markers.csv \
  -o primer_designs.csv \
  --include_fallback
```

---

## BLAST-assisted ranking

When BLAST screening is enabled:

- only the best `preblast_top_n` local candidates per marker are sent to BLAST
- BLAST uses `A1_Core`, `A2_Core`, and `Common_Core`
- final ranking combines local quality with off-target hit summaries

Without BLAST:

- output rows use `Ranking_mode = local_only`
- fallback-only outputs use `fallback_local_only`

With BLAST:

- output rows use `Ranking_mode = local_plus_blast`
- fallback-only outputs use `fallback_local_plus_blast`

### Example: keep 5 final rows per marker, but BLAST-screen the top 10 local rows first

```bash
python KASP_design.py \
  -g genome.fa \
  -m target_markers.csv \
  -o primer_designs.csv \
  --background_vcf input.vcf.gz \
  --blast_use_genome \
  --top_n 5 \
  --preblast_top_n 10
```

---

## Command-line options

## Required arguments

- `-g, --genome_path` — reference FASTA path
- `-m, --marker_csv` — marker CSV path
- `-o, --output_csv` — output CSV path

---

## Background masking options

- `--background_vcf` — optional VCF/VCF.gz of background variants
- `--background_csv` — optional CSV of background variants

---

## Local design settings

- `--allele_primer_min` — default `18`
- `--allele_primer_max` — default `30`
- `--common_primer_min` — default `18`
- `--common_primer_max` — default `30`
- `--tm_min` — default `57.0`
- `--tm_max` — default `65.0`
- `--tm_tolerance` — default `2.0`
- `--product_min` — default `60`
- `--product_max` — default `150`
- `--scan_window` — default `200`
- `--min_gc` — default `30.0`
- `--max_gc` — default `70.0`
- `--max_thermo_tm` — default `45.0`
- `--max_candidates_per_orientation` — default `200`

---

## Output retention settings

- `--top_n` — default `5`; maximum number of written rows per marker
- `--include_fallback` — write fallback rows when no `PASS` rows exist for that marker
- `--preblast_top_n` — default `10`; with BLAST enabled, only this many best local candidates per marker are BLAST-screened
- `--output_mode` — `concise` or `full`; default `concise`
- `--failed_output_csv` — optional path for the companion failed-marker report; if omitted, the script uses `<output stem>_failed.csv` when needed

---

## BLAST options

- `--blast_db` — existing BLAST database prefix or FASTA path
- `--blast_use_genome` — use `--genome_path` as the BLAST reference
- `--blast_db_dir` — optional directory for an auto-built BLAST database
- `--blastn_path` — default `blastn`
- `--makeblastdb_path` — default `makeblastdb`
- `--blast_task` — default `blastn-short`
- `--blast_word_size` — default `7`
- `--blast_max_target_seqs` — default `200`
- `--blast_num_threads` — default `1`

---

## Primer3 thermodynamic conditions

These parameters affect Primer3 thermodynamic calculations:

- `--mv_conc` — default `50.0` mM
- `--dv_conc` — default `1.5` mM
- `--dntp_conc` — default `0.6` mM
- `--dna_conc` — default `50.0` nM
- `--temp_c` — default `37.0` C
- `--max_loop` — default `30`

If your lab chemistry differs substantially from these defaults, consider adjusting them to better reflect your assay conditions.

---

## Output CSV structure

The output contains only **selected rows**. A single marker can contribute more than one row, up to `top_n`.

By default, `KASP_design.py` now writes a **concise** output schema. This keeps the design logic just as strict, but reduces the number of reported indicators after the main assay-definition and basic metric fields.

## Default concise output fields

- `Chr`
- `Position`
- `Ref`
- `Alt`
- `FASTA_base`
- `FASTA_allele_annotation`
- `TASSEL_major_allele_as_REF`
- `Orientation`
- `A1_Primer`
- `A2_Primer`
- `Common_Primer`
- `A1_Core`
- `A2_Core`
- `Common_Core`
- `Allele_specific_span_start`
- `Allele_specific_span_end`
- `Common_primer_span_start`
- `Common_primer_span_end`
- `A1_core_tm`
- `A2_core_tm`
- `Common_core_tm`
- `A1_core_gc`
- `A2_core_gc`
- `Common_core_gc`
- `Product_size`
- `passed`
- `worst_struct_tm`
- `max_tm_delta`
- `Local_score`
- `Ranking_mode`
- `Rank`
- `Design_status`

In other words, the default concise mode keeps all marker identity, FASTA check, primer sequence, primer span, and basic assay metric columns through `Product_size`, then condenses the downstream diagnostics to a small summary.

When BLAST screening is enabled, the concise output also adds:

- `Total_exact_full_length_off_target_hits`
- `Total_off_target_3prime_hits`
- `Worst_offtarget_3prime_bitscore`

Use `--output_mode full` if you want the full diagnostic output described below and in the appendix.

## Companion failed-marker report

Markers that do not receive any written assay rows are written to a separate CSV. By default, the path is derived from the main output path:

- main output: `primer_designs.csv`
- failed-marker report: `primer_designs_failed.csv`

The failed-marker report includes:

- `Chr`
- `Position`
- `Ref`
- `Alt`
- `FASTA_base`
- `FASTA_allele_annotation`
- `TASSEL_major_allele_as_REF`
- `Failure_reason`
- `Failure_reason_counts`
- `PASS_candidates_found`
- `FALLBACK_candidates_found`

`Failure_reason` is the dominant reason that no assay row was written for that marker, such as `missing_chromosome_in_fasta`, `position_out_of_range`, `ref_alt_do_not_match_fasta`, `no_candidate_found`, `thermo_fallback`, or the most common local primer-rejection category.

`Failure_reason_counts` is a semicolon-delimited summary of the internal local rejection counts when those counts exist, for example `forward_allele_specific_tm=12;reverse_allele_specific_overlaps_variant=8`.

## Marker and FASTA annotation fields

- `Chr`
- `Position`
- `Ref`
- `Alt`
- `FASTA_base`
- `FASTA_allele_annotation`
- `TASSEL_major_allele_as_REF`

## Assay definition fields

- `Orientation`
- `A1_Primer`
- `A2_Primer`
- `Common_Primer`
- `A1_Core`
- `A2_Core`
- `Common_Core`

## Genomic span fields

- `Allele_specific_span_start`
- `Allele_specific_span_end`
- `Common_primer_span_start`
- `Common_primer_span_end`

## Primer metric fields

- `A1_core_tm`
- `A2_core_tm`
- `Common_core_tm`
- `A1_core_gc`
- `A2_core_gc`
- `Common_core_gc`
- `Product_size`

## Structure metric fields

- `passed`
- `worst_struct_tm`
- `max_tm_delta`
- `A1_hairpin_tm`
- `A2_hairpin_tm`
- `Common_hairpin_tm`
- `A1_homodimer_tm`
- `A2_homodimer_tm`
- `Common_homodimer_tm`
- `A1_Common_heterodimer_tm`
- `A2_Common_heterodimer_tm`
- `A1_A2_heterodimer_tm`

## Ranking and status fields

- `Local_score`
- `Local_score_tuple`
- `Ranking_mode`
- `Rank`
- `Combined_rank_key`
- `Design_status`

## Additional BLAST fields

When BLAST screening is enabled, the output also includes:

- `A1_off_target_3prime_hits`
- `A2_off_target_3prime_hits`
- `Common_off_target_3prime_hits`
- `A1_exact_full_length_off_target_hits`
- `A2_exact_full_length_off_target_hits`
- `Common_exact_full_length_off_target_hits`
- `A1_best_offtarget_3prime_bitscore`
- `A2_best_offtarget_3prime_bitscore`
- `Common_best_offtarget_3prime_bitscore`
- `Total_off_target_3prime_hits`
- `Total_exact_full_length_off_target_hits`
- `Worst_offtarget_3prime_bitscore`

---

## How to read the FASTA annotation

### `REF matches FASTA`
The marker agrees with the reference genome at the target coordinate.

### `ALT matches FASTA; FASTA carries ALT allele`
The FASTA contains the ALT allele rather than the REF allele listed in the marker file. These markers are retained and annotated.

### `Neither REF nor ALT matches FASTA`
The marker is rejected. This usually indicates one of the following:

- wrong genome build
- shifted coordinates
- inconsistent chromosome naming
- marker file problem

---

## Terminal summary

At the end of each run, the script prints a summary including:

- `Output rows written`
- `Markers without written output rows`
- `Reference mismatches against FASTA`
- `ALT matches FASTA and likely treated as REF by TASSEL`
- `Outcome summary`
- the resolved BLAST database path when BLAST screening is enabled
- `BLAST database was auto-built for this run.` when applicable
- `Output saved to ...`
- `Failed-marker report saved to ...` when a companion report is written

### Common `Outcome summary` categories

- `designed`
- `thermo_fallback`
- `missing_chromosome_in_fasta`
- `position_out_of_range`
- `ref_alt_do_not_match_fasta`
- `*_tm`
- `*_gc`
- `*_overlaps_variant`
- `ambiguous_base_*`
- `no_candidate_found`

### Important interpretation note

`Reference mismatches against FASTA` is counted across processed markers, not only written output rows. In the current code, it increases whenever the marker is not annotated as `REF matches FASTA`, and markers that never reach allele-level FASTA annotation may also contribute.

---

## Practical examples

## Example 1: simplest possible run

```bash
python KASP_design.py \
  -g genome.fa \
  -m markers.csv \
  -o kasp_out.csv
```

Use this when you already trust the marker file and do not need background masking or BLAST.

## Example 2: recommended design run with variant masking

```bash
python KASP_design.py \
  -g genome.fa \
  -m markers.csv \
  -o kasp_out.csv \
  --background_vcf all_variants.vcf.gz  #(or --backgroud_csv)
```

Use this when your markers were selected from a larger variant call set and you want to avoid adjacent known variants inside primer-binding regions.

## Example 3: stricter amplicon and structure filtering

```bash
python KASP_design.py \
  -g genome.fa \
  -m markers.csv \
  -o kasp_out.csv \
  --background_vcf all_variants.vcf.gz \   #(or --backgroud_csv)
  --product_min 70 \
  --product_max 110 \
  --max_thermo_tm 38
```

Use this when you want to be more conservative and accept fewer, cleaner designs.

## Example 4: keep more candidate rows per marker

```bash
python KASP_design.py \
  -g genome.fa \
  -m markers.csv \
  -o kasp_out.csv \
  --background_vcf all_variants.vcf.gz \         #(or --backgroud_csv)
  --top_n 10
```

Use this when you want to manually inspect several ranked assays per SNP.

## Example 5: include fallback rows

```bash
python KASP_design.py \
  -g genome.fa \
  -m markers.csv \
  -o kasp_out.csv \
  --background_vcf all_variants.vcf.gz \         #(or --backgroud_csv)
  --include_fallback
```

Use this when you prefer to see weak but still usable backup designs rather than dropping those markers entirely.

## Example 6: BLAST-assisted ranking with genome-based auto database creation

```bash
python KASP_design.py \
  -g genome.fa \
  -m markers.csv \
  -o kasp_out.csv \
  --background_vcf all_variants.vcf.gz \     #(or --backgroud_csv)
  --blast_use_genome \ 
  --top_n 5 \
  --product_min 60 \
  --product_max 120 \
  --max_thermo_tm 40
  --preblast_top_n 10 \
  --blast_num_threads 4
```

Use this when you want local candidate ranking plus off-target visibility across the genome.

---

## Recommended workflow

1. Convert your VCF to a clean marker CSV with `vcf_to_kasp_csv.py`.
2. Validate marker positions and REF alleles against the reference with `check_markers_against_fasta.py`.
3. Run `KASP_design.py` using the original VCF as `--background_vcf` when available.
4. Keep multiple ranked assays per marker if you want choice for wet-lab testing.
5. Add BLAST screening when specificity is a concern or when you want BLAST-assisted ranking.

---

## Troubleshooting

## Output is empty or unexpectedly small

Check these first:

- FASTA chromosome names really match the marker table after `chr` normalization
- marker positions are 1-based
- the marker file contains only simple SNPs, not indels or multiallelic rows
- at least one of `REF` or `ALT` matches the FASTA
- the background VCF or CSV is not so dense that nearly all primer candidates overlap another variant
- `--include_fallback` is off by default

## `Output rows written: 0`

This means no marker produced a selected row under the current settings. The most useful next step is to inspect the terminal `Outcome summary`, because the script records the dominant failure reason for each marker.

## Many markers fail as `ref_alt_do_not_match_fasta`

Possible causes:

- wrong reference genome version
- positions are not 1-based
- chromosome names do not match after normalization
- marker CSV came from a different coordinate system

Run `check_markers_against_fasta.py` first to confirm.

## Many markers fail as `*_overlaps_variant`

This usually means the surrounding region is highly polymorphic under your strict masking rule. You may need to:

- choose a different SNP set
- reduce background density if it contains unwanted low-confidence variants
- increase the number of candidate regions scanned only if the code supports broader search

## BLAST errors

Make sure:

- `blastn` is installed and visible in `PATH`
- `makeblastdb` is installed if auto-building is enabled
- the FASTA is readable
- the BLAST database path points to the database prefix, not only to one database sidecar file

---

## Design philosophy

This repository is intentionally conservative.

It prefers:

- rejecting primers that overlap known nearby variants
- retaining multiple ranked designs per SNP
- using Primer3 for thermodynamic checks
- optionally adding BLAST-based off-target visibility

rather than simply returning the first sequence pair that appears to work.

That makes the output more selective, but generally better suited for downstream assay testing.



## a real case
```
python /home/cflzxc/software/Batch_KASP/KASP_design.py \
  -g /output/genomic/plant/Pisum/sativum/Genome/Pisum_sativum_v1a.fa \
  -m markers_for_kasp.csv \
  -o primer_designs.csv \
  --background_csv ALL_SNPs.csv
  --blast_use_genome \ 
  --top_n 5 \
  --product_min 60 \
  --product_max 120 \
  --max_thermo_tm 40
  --preblast_top_n 10 \
  --blast_num_threads 4


  one_line 
  python /home/cflzxc/software/Batch_KASP/KASP_design.py -g /output/genomic/plant/Pisum/sativum/Genome/Pisum_sativum_v1a.fa -m markers_for_kasp.csv -o primer_designs.csv --background_csv ALL_SNPs.csv --blast_use_genome --top_n 5 --product_min 60 --product_max 120 --max_thermo_tm 40 --preblast_top_n 10 --blast_num_threads 4 --include_fallback
```

---

## Appendix: full output columns from `KASP_design.py`

The list below annotates the columns written when `KASP_design.py` is run with `--output_mode full`. Some columns appear only when BLAST screening is enabled.

### General notes

- `A1` is the REF-linked assay primer in the current code
- `A2` is the ALT-linked assay primer in the current code
- all genomic coordinates are 1-based
- all `*_tm` values are Primer3-derived temperatures in degrees C
- all `*_gc` values are GC percentages
- `Common_Primer` and `Common_Core` are the same sequence in the current implementation

### Marker identity, FASTA checks, and primer sequences

| Column | Meaning |
|---|---|
| `Chr` | Chromosome or contig name from the normalized marker table. |
| `Position` | Target SNP position in 1-based genomic coordinates. |
| `Ref` | REF allele from the marker CSV after normalization to uppercase. |
| `Alt` | ALT allele from the marker CSV after normalization to uppercase. |
| `FASTA_base` | Base found in the reference FASTA at `Chr:Position`. |
| `FASTA_allele_annotation` | FASTA consistency label. Current values are `REF matches FASTA`, `ALT matches FASTA; FASTA carries ALT allele`, or `Neither REF nor ALT matches FASTA`. |
| `TASSEL_major_allele_as_REF` | Current TASSEL-style interpretation flag from the script. `No` when REF matches FASTA, `Yes` when ALT matches FASTA, and `Unknown` when the marker could not be cleanly aligned to the FASTA. |
| `Orientation` | Assay layout chosen for that row. Current values are `allele_specific_forward_common_reverse` or `allele_specific_reverse_common_forward`. |
| `A1_Primer` | Full REF-specific allele primer sequence, including the FAM tail plus `A1_Core`. |
| `A2_Primer` | Full ALT-specific allele primer sequence, including the VIC tail plus `A2_Core`. |
| `Common_Primer` | Full common primer sequence used with the allele-specific pair. In the current code this is the untailored common primer core. |
| `A1_Core` | REF-specific allele primer core sequence, without the FAM tail. |
| `A2_Core` | ALT-specific allele primer core sequence, without the VIC tail. |
| `Common_Core` | Common primer core sequence. In the current code this is identical to `Common_Primer`. |

### Primer genomic spans and basic assay metrics

| Column | Meaning |
|---|---|
| `Allele_specific_span_start` | Start coordinate of the allele-specific primer binding span. |
| `Allele_specific_span_end` | End coordinate of the allele-specific primer binding span. This span includes the target SNP. |
| `Common_primer_span_start` | Start coordinate of the common primer binding span. |
| `Common_primer_span_end` | End coordinate of the common primer binding span. |
| `A1_core_tm` | Melting temperature of `A1_Core`. |
| `A2_core_tm` | Melting temperature of `A2_Core`. |
| `Common_core_tm` | Melting temperature of `Common_Core`. |
| `A1_core_gc` | GC percentage of `A1_Core`. |
| `A2_core_gc` | GC percentage of `A2_Core`. |
| `Common_core_gc` | GC percentage of `Common_Core`. |
| `Product_size` | Expected amplicon size in base pairs for the primer pair. |

### Thermodynamic and structure-check columns

| Column | Meaning |
|---|---|
| `passed` | Raw pass/fail flag before final row selection. It is `True` only when `Product_size` is inside the configured range and `worst_struct_tm <= max_thermo_tm`. |
| `worst_struct_tm` | Worst, meaning highest, structure Tm observed across all checked hairpins, homodimers, and heterodimers for the candidate set. |
| `max_tm_delta` | Largest absolute Tm difference among `A1_Core`, `A2_Core`, and `Common_Core`. |
| `A1_hairpin_tm` | Hairpin Tm for `A1_Primer`. |
| `A2_hairpin_tm` | Hairpin Tm for `A2_Primer`. |
| `Common_hairpin_tm` | Hairpin Tm for `Common_Primer`. |
| `A1_homodimer_tm` | Homodimer Tm for `A1_Primer`. |
| `A2_homodimer_tm` | Homodimer Tm for `A2_Primer`. |
| `Common_homodimer_tm` | Homodimer Tm for `Common_Primer`. |
| `A1_Common_heterodimer_tm` | Heterodimer Tm between `A1_Primer` and `Common_Primer`. |
| `A2_Common_heterodimer_tm` | Heterodimer Tm between `A2_Primer` and `Common_Primer`. |
| `A1_A2_heterodimer_tm` | Heterodimer Tm between `A1_Primer` and `A2_Primer`. |

### Local ranking and final selection columns

| Column | Meaning |
|---|---|
| `Local_score` | Scalar local ranking score computed from structure Tm, Tm balance, product-size deviation, GC penalty, and Tm-center penalty. Lower is better. |
| `Local_score_tuple` | Pipe-delimited tuple used for the script's local candidate ordering: `worst_struct_tm|max_tm_delta|product_size_distance_from_midpoint|gc_penalty|tm_center_penalty`. Lower tuples rank better. |
| `Ranking_mode` | How the final row was ranked before writing. Current values are `local_only`, `local_plus_blast`, `fallback_local_only`, and `fallback_local_plus_blast`. |
| `Rank` | 1-based rank among the rows written for that SNP in the final output. |
| `Combined_rank_key` | The final sort key written as text. Without BLAST it is the string form of `Local_score`. With BLAST it is the pipe-joined tuple `Total_exact_full_length_off_target_hits|Total_off_target_3prime_hits|-Worst_offtarget_3prime_bitscore|Local_score|worst_struct_tm|max_tm_delta`. |
| `Design_status` | Final written status for the row: `PASS` or `FALLBACK`. `FALLBACK` rows are only written when `--include_fallback` is enabled and the marker has no `PASS` candidates. |

### BLAST-only off-target columns

These columns are added only when BLAST screening is enabled with `--blast_db` or `--blast_use_genome`.

| Column | Meaning |
|---|---|
| `A1_off_target_3prime_hits` | Number of off-target BLAST hits for `A1_Core` that cover the primer 3' end and are near full length under the script's current filters. |
| `A2_off_target_3prime_hits` | Number of off-target BLAST hits for `A2_Core` that cover the primer 3' end and are near full length under the script's current filters. |
| `Common_off_target_3prime_hits` | Number of off-target BLAST hits for `Common_Core` that cover the primer 3' end and are near full length under the script's current filters. |
| `A1_exact_full_length_off_target_hits` | Number of exact full-length off-target BLAST matches for `A1_Core`. |
| `A2_exact_full_length_off_target_hits` | Number of exact full-length off-target BLAST matches for `A2_Core`. |
| `Common_exact_full_length_off_target_hits` | Number of exact full-length off-target BLAST matches for `Common_Core`. |
| `A1_best_offtarget_3prime_bitscore` | Best, meaning highest, BLAST bitscore among off-target 3'-covering near-full-length hits for `A1_Core`. |
| `A2_best_offtarget_3prime_bitscore` | Best, meaning highest, BLAST bitscore among off-target 3'-covering near-full-length hits for `A2_Core`. |
| `Common_best_offtarget_3prime_bitscore` | Best, meaning highest, BLAST bitscore among off-target 3'-covering near-full-length hits for `Common_Core`. |
| `Total_off_target_3prime_hits` | Sum of `A1_off_target_3prime_hits`, `A2_off_target_3prime_hits`, and `Common_off_target_3prime_hits`. |
| `Total_exact_full_length_off_target_hits` | Sum of `A1_exact_full_length_off_target_hits`, `A2_exact_full_length_off_target_hits`, and `Common_exact_full_length_off_target_hits`. |
| `Worst_offtarget_3prime_bitscore` | Maximum of the three per-primer best off-target 3' bitscores. In practice this is the strongest off-target 3'-covering bitscore seen among the three primers. |
