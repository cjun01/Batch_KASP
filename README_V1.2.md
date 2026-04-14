# Batch_KASP

`Batch_KASP` is a batch KASP assay design workflow for **biallelic SNP markers and simple biallelic indels**. It takes a **reference genome FASTA** plus a **marker CSV/TSV** and returns one or more ranked KASP assay designs per target marker.

This README is updated for the indel-capable script:

- `KASP_design_V1.2.py`

If you later rename that file back to `KASP_design.py`, you can keep this README and just update the script name in the examples.

The main design engine:

- reads a reference FASTA or FASTA.gz
- validates and normalizes marker input
- supports SNPs plus simple insertions/deletions
- accepts a human-readable deletion input style such as `13969-13971,GAT,-`
- supports both common KASP assay orientations
- rejects primers that overlap nearby known variants by default
- can optionally tolerate a limited shared non-target variant overlap in the common region of both A1 and A2, with explicit warning columns in the output
- uses **Primer3 thermodynamic checks** for hairpins and dimer interactions
- keeps multiple valid candidates per marker
- ranks candidates locally
- can optionally re-rank top candidates using **BLAST off-target screening**

The repository may also include helper utilities such as:

- `vcf_to_kasp_csv.py` — convert a VCF/VCF.gz into a marker CSV
- `check_markers_against_fasta.py` — verify marker coordinates and alleles against the reference FASTA

---

## Repository contents

```text
Batch_KASP/
├── KASP_design_V1.2.py
├── vcf_to_kasp_csv.py
├── check_markers_against_fasta.py
├── requirements.txt
└── README.md
```

---

## What the pipeline does

For each marker, `KASP_design_V1.2.py` performs the following steps:

1. **Loads the reference genome**
   - accepts FASTA or FASTA.gz
   - normalizes chromosome names so `chr1` and `1` are treated as the same chromosome

2. **Reads and cleans the marker table**
   - requires `Chr`, `position`, `ref`, and `alt`
   - treats column names case-insensitively
   - accepts **CSV/TSV marker tables only** for `--marker_csv`
   - accepts supported SNP and indel formats inside that table
   - can split comma-separated ALT alleles in the `alt` column into separate biallelic design targets
   - removes exact duplicate marker rows after normalization
   - rejects raw `.vcf` and `.vcf.gz` files as `--marker_csv` input

3. **Checks allele consistency against the FASTA**
   - `REF matches FASTA`
   - `ALT matches FASTA; FASTA carries ALT allele`
   - `REF does not match FASTA at marker span`

   Notes:
   - SNPs can be retained when the FASTA carries the ALT allele instead of the REF allele.
   - Indels are currently validated primarily by requiring the **REF allele span** to match the FASTA at the marker span.

4. **Builds a strict variant mask for primer placement**
   - all target markers are masked automatically
   - optional background variants from a VCF or CSV are also masked
   - the current target marker span is still allowed inside its own allele-specific primer region
   - by default, any other overlapping known variant causes that primer candidate to be rejected
   - optionally, the script can tolerate a very limited exception where the **same non-target variant position(s)** overlap both `A1` and `A2` in their shared region and remain sufficiently far from both 3' ends; accepted rows are flagged in the output

5. **Searches supported KASP assay layouts**
   - for SNPs:
     - **Orientation 1**: allele-specific forward primers + common reverse primer
     - **Orientation 2**: allele-specific reverse primers + common forward primer
   - for indels:
     - **Orientation 1**: indel allele-specific forward primers + common reverse primer
     - **Orientation 2**: indel allele-specific reverse primers + common forward primer

6. **Adds the standard reporter tails**
   - `A1_Primer = FAM tail + A1_Core`
   - `A2_Primer = VIC tail + A2_Core`

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
   - multiple ranked candidates can be retained per marker

9. **Optionally performs BLAST off-target screening**
   - BLAST is run on `A1_Core`, `A2_Core`, and `Common_Core`
   - top local candidates are re-ranked using off-target hit summaries
   - indel-aware per-primer spans are used when BLAST queries are built

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

### 1) Reference FASTA

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

### 2) Marker CSV / TSV

Argument:

```bash
-m / --marker_csv
```

Required columns:
- `Chr`
- `position`
- `ref`
- `alt`

Column matching is case-insensitive. The marker file may be comma-separated or tab-separated, but it must remain a **marker table**, not a raw VCF file.

### Supported marker formats

#### A. SNP

```csv
Chr,position,ref,alt
1LG6,13969,G,T
```

Rules:
- `position` is a single **1-based** coordinate
- `ref` and `alt` are single A/C/G/T bases

#### B. Anchored insertion

```csv
Chr,position,ref,alt
1LG6,1300,G,GATC
```

This means the reference allele is `G` and the alternate allele is `GATC` at the anchor position.

#### C. Anchored deletion

```csv
Chr,position,ref,alt
1LG6,1300,GATC,G
```

This means the reference allele is `GATC` and the alternate allele is `G` at the anchor position.

#### C2. Comma-separated ALT alleles in the marker table

```csv
Chr,position,ref,alt
ChrA01,30423823,TAA,"TCTTAAA,TCTTAAAA"
```

or, if saved as a TSV, the same logical columns with the `alt` cell containing:

```text
TCTTAAA,TCTTAAAA
```

The script splits this one marker-table row into separate biallelic design targets internally:
- `TAA -> TCTTAAA`
- `TAA -> TCTTAAAA`

This is useful when you want to keep the older marker CSV/TSV workflow but still preserve multiple ALT alleles from upstream processing.

#### D. Human-readable deletion range

```csv
Chr,position,ref,alt
1LG6,13969-13971,GAT,-
```

This means the reference genome contains `GAT` from positions `13969` to `13971`, and the alternate allele deletes that block.

For this human-readable deletion style:
- `position` must be a range like `start-end`
- `ref` must match the deleted block length
- `alt` can be written as:
  - empty
  - `-`
  - `.`
  - `DEL`
  - `DELETION`

Internally, the script keeps:
- `Position` = the start coordinate of the deleted block
- `Input_position` = the original token, such as `13969-13971`

### Currently unsupported marker formats

These are rejected:
- equal-length multibase substitutions such as `AT -> GC`
- range rows that are not deletion-style
- numeric rows with an empty ALT allele
- raw `.vcf` or `.vcf.gz` files passed directly as `--marker_csv`

### Notes on indel input style

- For **human readability**, deletion ranges like `13969-13971,GAT,-` are convenient.
- For **insertions**, the script currently expects **anchored input**, not a human-readable range style.
- If your indels come from a VCF, convert or copy them into anchored CSV/TSV rows. The design script no longer accepts a raw VCF file directly as `--marker_csv`.

---

### 3) Optional background variants

You can provide known nearby variants that should be avoided in primer-binding regions.

**Background VCF:**

```bash
--background_vcf input.vcf.gz
```

Accepted: `.vcf`, `.vcf.gz`. The masked span uses the length of the VCF `REF` allele.

**Background CSV:**

```bash
--background_csv background.csv
```

Must contain at least `Chr` and `position`.

Behavior:
- if a `ref` column is present, the masked span uses its length
- otherwise a single base is masked
- deletion-style ranges such as `13969-13971` are also accepted in the background CSV and will mask the full span

**Important masking rule**

By default, the script is intentionally strict:
- the current target marker span is allowed for its own assay
- any **other** masked variant overlapping a primer binding region causes rejection

Optional exception:
- if you enable `--allow_shared_allele_primer_variants`, the script may keep a candidate when the **same non-target position(s)** overlap both `A1` and `A2`, the overlap count is within the configured limit, and the overlap remains far enough from both 3' ends
- accepted rows are flagged with warning columns such as `Shared_variant_tolerance_used`, `Shared_variant_positions`, and `Shared_variant_warning`

This means the pipeline still prefers **cleaner assays** over maximizing the number of returned designs, but now allows a controlled and transparent compromise for difficult loci.

---

## Workflow overview

A recommended workflow is:

1. Convert a VCF to a clean marker CSV
2. Check that marker coordinates and alleles match the FASTA
3. Run local KASP design
4. Optionally add BLAST off-target screening

---

## Step 1: Prepare the marker CSV

You can prepare a marker CSV manually or convert from VCF.

### Convert from VCF with `vcf_to_kasp_csv.py`

Use the helper script to convert a `.vcf` or `.vcf.gz` file into the marker CSV format required by the design script. This helper section is retained intentionally because it is still useful for SNP workflows:

```bash
python vcf_to_kasp_csv.py \
  -i input.vcf.gz \
  -o target_markers.csv
```

Keep only markers for a specific sample:

```bash
python vcf_to_kasp_csv.py \
  -i input.vcf.gz \
  -o target_markers.csv \
  -s SampleName
```

Keep only non-reference genotypes (`0/1` and `1/1`) for that sample:

```bash
python vcf_to_kasp_csv.py \
  -i input.vcf.gz \
  -o target_markers.csv \
  -s SampleName \
  --genotype nonref
```

Skip contigs whose names contain a substring such as `unitig`:

```bash
python vcf_to_kasp_csv.py \
  -i input.vcf.gz \
  -o target_markers.csv \
  --exclude-contig-substring unitig
```

Show all options:

```bash
python vcf_to_kasp_csv.py --help
```

Notes:
- the output columns are `Chr,position,ref,alt`
- the helper writes only simple biallelic SNPs
- indels and most multi-allelic sites are skipped by the helper and still need to be added manually or by extending the converter
- if you provide `-s/--sample` and omit `--genotype`, the default genotype filter is `alt` (`1/1`)

If you need indel rows, you may need to add them manually or extend the converter.

Examples of accepted rows:

```csv
Chr,position,ref,alt
1LG6,13969,G,T
1LG6,1300,G,GATC
1LG6,1300,GATC,G
1LG6,13969-13971,GAT,-
ChrA01,30423823,TAA,"TCTTAAA,TCTTAAAA"
```

---

## Step 2: Check markers against the FASTA

Before design, verify that the marker coordinates and REF alleles are consistent with the reference.

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
- cases where `ALT`, not `REF`, matches the FASTA for SNPs
- indel rows whose REF span does not actually match the FASTA

You should inspect the resulting report before running large design jobs.

---

## Step 3: Run local KASP design

Minimal local design example:

```bash
python KASP_design_V1.2.py \
  -g genome.fa \
  -m target_markers.csv \
  -o primer_designs.csv
```

By default, the script writes a **concise output table**: it keeps the internal screening and ranking logic, but only writes a smaller set of decision-useful columns. Use `--output_mode full` if you want the full diagnostic table.

`--marker_csv` must still point to a marker CSV/TSV, not a raw VCF file. Use `vcf_to_kasp_csv.py` for SNP extraction or prepare indel rows manually in the same table format.

When shared-variant tolerance is enabled and used for a returned row, the concise output also reports warning fields describing the tolerated overlap.

If a marker does not receive any written assay row, the script also writes a companion failed-marker report by default using the output stem plus `_failed.csv`. For example, `primer_designs.csv` is paired with `primer_designs_failed.csv`.

### Example: include a background VCF

```bash
python KASP_design_V1.2.py \
  -g genome.fa \
  -m target_markers.csv \
  -o primer_designs.csv \
  --background_vcf input.vcf.gz
```

### Example: include a background CSV

```bash
python KASP_design_V1.2.py \
  -g genome.fa \
  -m target_markers.csv \
  -o primer_designs.csv \
  --background_csv background.csv
```

### Example: stricter local filters

```bash
python KASP_design_V1.2.py \
  -g genome.fa \
  -m target_markers.csv \
  -o primer_designs.csv \
  --background_vcf input.vcf.gz \
  --product_min 60 \
  --product_max 120 \
  --max_thermo_tm 40
```

### Example: allow a controlled shared-variant exception in A1/A2

Use this only for difficult loci where the same unavoidable non-target SNP position overlaps the shared region of both allele-specific primers.

```bash
python KASP_design_V1.2.py \
  -g genome.fa \
  -m target_markers.csv \
  -o primer_designs.csv \
  --background_vcf input.vcf.gz \
  --allow_shared_allele_primer_variants \
  --max_shared_allele_variant_positions 1 \
  --min_shared_allele_variant_distance_3prime 5
```

This mode still rejects most overlap situations. It only tolerates shared non-target positions that:
- hit both `A1` and `A2`
- stay within the configured maximum count
- remain at least the configured distance away from both 3' ends

Accepted rows are clearly flagged in the output.

---

## Step 4: Add BLAST off-target screening

BLAST screening is optional and is mainly useful when the genome is large or repetitive, duplicated loci are a concern, or you want BLAST-assisted ranking among already acceptable local candidates.

### Option A: use an existing BLAST database

```bash
python KASP_design_V1.2.py \
  -g genome.fa \
  -m target_markers.csv \
  -o primer_designs.csv \
  --background_vcf input.vcf.gz \
  --blast_db /path/to/blast_db_prefix
```

### Option B: use the genome FASTA directly and auto-build the BLAST database

```bash
python KASP_design_V1.2.py \
  -g genome.fa \
  -m target_markers.csv \
  -o primer_designs.csv \
  --background_vcf input.vcf.gz \
  --blast_use_genome
```

---

## Candidate generation and ranking

The script does **not** stop at the first acceptable assay. Instead, it can keep multiple acceptable candidates per marker and rank them.

### Local candidate generation

For each marker:
- supported orientations are explored
- allele-specific and common primers must satisfy current length, Tm, and GC constraints
- by default, any primer overlapping a masked background variant is rejected, except for the marker's own target span
- when `--allow_shared_allele_primer_variants` is enabled, a limited shared-overlap exception can be applied only when the same non-target position(s) overlap both `A1` and `A2` and remain sufficiently far from both 3' ends
- Primer3 structure checks are computed for each surviving candidate pair
- each candidate gets a local ranking record

### SNP logic

For SNPs:
- the two allele-specific primers differ at the allele-discriminating end
- the common primer is chosen on the opposite shared flank

### Indel logic

For simple indels:
- `A1` is the **REF-linked** allele-specific primer
- `A2` is the **ALT-linked** allele-specific primer
- the design engine builds REF and ALT local haplotypes around the target
- allele-specific primers are generated to discriminate the REF and ALT junction structures
- the common primer is chosen from the opposite shared flank

For human-readable deletion rows like:

```csv
1LG6,13969-13971,GAT,-
```

this corresponds conceptually to:
- `A1` binding a sequence that includes the deleted REF block
- `A2` binding the ALT haplotype that skips that block and spans the new junction

### Local ranking logic

Candidates are locally ordered by `Local_score_tuple`, which stores:
1. `worst_struct_tm`
2. `max_tm_delta`
3. distance from the ideal product-size midpoint
4. GC penalty relative to 50%
5. Tm penalty relative to the midpoint of `tm_min` and `tm_max`

The script writes a scalar `Local_score`, but the primary local ordering is driven by this tuple.

When shared-variant tolerance is used, the script adds an extra penalty so cleaner assays still rank ahead of tolerated-warning designs when otherwise comparable.

---

## `PASS` vs `FALLBACK`

Each candidate is first assessed locally. The boolean field `passed` is `True` only if both conditions are satisfied:
- `product_min <= Product_size <= product_max`
- `worst_struct_tm <= max_thermo_tm`

Selected output rows use `Design_status`:
- `PASS` — selected from candidates that satisfy the local pass criteria
- `FALLBACK` — selected from weaker candidates only when no `PASS` candidates exist

A row can still be `PASS` even when `Shared_variant_tolerance_used = Yes`, as long as it meets the thermodynamic and product-size thresholds. The warning fields tell you that a controlled overlap exception was used.

By default, only `PASS` rows are written. To include fallback rows when a marker has no passing candidate:

```bash
python KASP_design_V1.2.py \
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

---

## Command-line options

**Required arguments**
- `-g, --genome_path` — reference FASTA path
- `-m, --marker_csv` — marker CSV path
- `-o, --output_csv` — output CSV path

**Background masking options**
- `--background_vcf` — optional VCF/VCF.gz of background variants
- `--background_csv` — optional CSV of background variants

**Local design settings**
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

**Shared-variant tolerance settings**
- `--allow_shared_allele_primer_variants` — enable a limited exception where the same non-target variant position(s) may overlap both `A1` and `A2` in their shared region, provided the overlap stays sufficiently far from both 3' ends; accepted rows are flagged
- `--max_shared_allele_variant_positions` — default `1`; maximum number of shared non-target variant positions allowed in both `A1` and `A2` when the tolerance flag is enabled
- `--min_shared_allele_variant_distance_3prime` — default `5`; minimum required distance in nucleotides from the 3' end of both `A1` and `A2` for any tolerated shared non-target variant

**Output retention settings**
- `--top_n` — default `5`; maximum number of written rows per marker
- `--min_common_primer_gap` — default `5`; require output rows for the same marker to have common-primer binding sites separated by at least this many bases
- `--include_fallback` — write fallback rows when no `PASS` rows exist for that marker
- `--preblast_top_n` — default `10`; with BLAST enabled, only this many best local candidates per marker are BLAST-screened
- `--output_mode` — `concise` or `full`; default `concise`
- `--failed_output_csv` — optional path for the companion failed-marker report
- `--tall_output_csv` — optional tall-format output with one primer per row

**BLAST options**
- `--blast_db` — existing BLAST database prefix or FASTA path
- `--blast_use_genome` — use `--genome_path` as the BLAST reference
- `--blast_db_dir` — optional directory for an auto-built BLAST database
- `--blastn_path` — default `blastn`
- `--makeblastdb_path` — default `makeblastdb`
- `--blast_task` — default `blastn-short`
- `--blast_word_size` — default `7`
- `--blast_max_target_seqs` — default `200`
- `--blast_num_threads` — default `1`

**Primer3 thermodynamic conditions**
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

By default, the script writes a **concise** output schema. Use `--output_mode full` if you want the full diagnostic table.

### Important output identity columns

- **`Chr`**: chromosome or contig name
- **`Position`**: normalized 1-based marker coordinate used internally
  - for deletion ranges, this is the **start** of the deleted block
- **`Input_position`**: original input token from the marker CSV
  - examples: `13969`, `1300`, `13969-13971`
- **`Ref`**: normalized REF allele
- **`Alt`**: normalized ALT allele

### FASTA check columns

- **`FASTA_base`**: for SNPs, this is the FASTA base at the marker; for indels, it reflects the FASTA sequence slice reported by the current code
- **`FASTA_allele_annotation`**: marker-to-FASTA consistency message
- **`TASSEL_major_allele_as_REF`**: `Yes`, `No`, or `Unknown`

### Primer sequence columns

- **`Orientation`**
  - SNP:
    - `allele_specific_forward_common_reverse`
    - `allele_specific_reverse_common_forward`
  - indel:
    - `indel_allele_specific_forward_common_reverse`
    - `indel_allele_specific_reverse_common_forward`
- **`A1_Primer`**: full primer to order for the REF-linked allele-specific primer
- **`A2_Primer`**: full primer to order for the ALT-linked allele-specific primer
- **`Common_Primer`**: full common primer to order
- **`A1_Core` / `A2_Core` / `Common_Core`**: genomic binding sequences without reporter tails

### Primer span columns

For SNPs, `A1` and `A2` usually occupy the same genomic span. For indels, they can occupy different effective spans.

- **`A1_span_start` / `A1_span_end`**: genomic span used for `A1_Core`
- **`A2_span_start` / `A2_span_end`**: genomic span used for `A2_Core`
- **`Allele_specific_span_start` / `Allele_specific_span_end`**: the union span across the two allele-specific primers
- **`Common_primer_span_start` / `Common_primer_span_end`**: genomic span for the common primer

### Basic assay metrics

- **`A1_core_tm` / `A2_core_tm` / `Common_core_tm`**
- **`A1_core_gc` / `A2_core_gc` / `Common_core_gc`**
- **`Product_size`**

For indels, `Product_size` is reported conservatively from the larger of the REF-linked and ALT-linked product calculations used internally.

### Scoring and ranking columns

- **`passed`**
- **`Exclusion_reason`**
- **`worst_struct_tm`**
- **`max_tm_delta`**
- **`Local_score`**
- **`Ranking_mode`**
- **`Rank`**
- **`Design_status`**

### BLAST columns

These appear only when BLAST is enabled:
- **`Total_exact_full_length_off_target_hits`**
- **`Total_off_target_3prime_hits`**
- **`Worst_offtarget_3prime_bitscore`**

### Tall-format output

If you use:

```bash
--tall_output_csv primer_sequences_tall.csv
```

then the script also writes a tall file with:
- `Primer_ID`
- `Primer_Sequence`

The `Primer_ID` format is:

```text
Chr_Position_rank<Rank>_<PrimerColumn>
```

---

## Companion failed-marker report

Markers that do not receive any written assay rows are written to a separate CSV. By default, the path is derived from the main output path:

- main output: `primer_designs.csv`
- failed-marker report: `primer_designs_failed.csv`

The failed-marker report includes:

- **`Chr`, `Position`, `Input_position`, `Ref`, `Alt`**: normalized and original marker details
- **`FASTA_base`, `FASTA_allele_annotation`, `TASSEL_major_allele_as_REF`**: FASTA alignment checks
- **`Failure_reason`**: primary reason no assay row was written
- **`Failure_reason_counts`**: detailed tally of internal candidate rejection reasons
- **`PASS_candidates_found` / `FALLBACK_candidates_found`**: counts of retained candidates before final writing

Examples of failure reasons include:
- `missing_chromosome_in_fasta`
- `position_out_of_range`
- `ref_alt_do_not_match_fasta`
- `equal_length_multibase_not_supported`
- `indel_forward_allele_specific_candidates`
- `indel_reverse_allele_specific_candidates`
- `no_candidate_found`

---

## How to read the FASTA annotation

### `REF matches FASTA`
The marker agrees with the reference genome at the target coordinate or target span.

### `ALT matches FASTA; FASTA carries ALT allele`
The FASTA contains the ALT allele rather than the REF allele listed in the marker file. This is mainly relevant to SNP-style markers.

### `REF does not match FASTA at marker span`
The marker is rejected. This usually indicates one of the following:
- wrong genome build
- shifted coordinates
- inconsistent chromosome naming
- marker file problem
- incorrectly normalized indel input

---

## Terminal summary

At the end of each run, the script prints a summary including:
- `Output rows written`
- `Output mode`
- `Markers without written output rows`
- `Reference mismatches against FASTA`
- `ALT matches FASTA and likely treated as REF by TASSEL`
- `Outcome summary`
- the resolved BLAST database path when BLAST screening is enabled
- `BLAST database was auto-built for this run.` when applicable
- `Output saved to ...`
- `Tall-format output saved to ...` when applicable
- `Failed-marker report saved to ...` when a companion report is written

---

## Practical examples

### Example 1: simplest possible run

```bash
python KASP_design_V1.2.py -g genome.fa -m markers.csv -o kasp_out.csv
```

### Example 2: design run with background masking

```bash
python KASP_design_V1.2.py \
  -g genome.fa \
  -m markers.csv \
  -o kasp_out.csv \
  --background_vcf all_variants.vcf.gz
```

### Example 3: include human-readable deletion rows

Marker CSV:

```csv
Chr,position,ref,alt
1LG6,13969,G,T
1LG6,13969-13971,GAT,-
```

Run:

```bash
python KASP_design_V1.2.py \
  -g genome.fa \
  -m markers.csv \
  -o kasp_out.csv
```

### Example 4: keep more candidate rows per marker

```bash
python KASP_design_V1.2.py \
  -g genome.fa \
  -m markers.csv \
  -o kasp_out.csv \
  --top_n 10
```

### Example 5: include fallback rows

```bash
python KASP_design_V1.2.py \
  -g genome.fa \
  -m markers.csv \
  -o kasp_out.csv \
  --include_fallback
```

### Example 6: BLAST-assisted ranking with auto database creation

```bash
python KASP_design_V1.2.py \
  -g genome.fa \
  -m markers.csv \
  -o kasp_out.csv \
  --background_vcf all_variants.vcf.gz \
  --blast_use_genome \
  --top_n 5 \
  --product_min 60 \
  --product_max 120 \
  --max_thermo_tm 40 \
  --preblast_top_n 10 \
  --blast_num_threads 4
```

### Example 7: allow a flagged shared-variant compromise design

```bash
python KASP_design_V1.2.py \
  -g genome.fa \
  -m markers.csv \
  -o kasp_out.csv \
  --background_vcf all_variants.vcf.gz \
  --allow_shared_allele_primer_variants \
  --max_shared_allele_variant_positions 1 \
  --min_shared_allele_variant_distance_3prime 5
```

### Example 8: write a tall primer list

```bash
python KASP_design_V1.2.py \
  -g genome.fa \
  -m markers.csv \
  -o kasp_out.csv \
  --tall_output_csv kasp_out_tall.csv
```

---

## Recommended workflow

1. Prepare a clean marker CSV/TSV containing supported SNP and indel rows.
2. Validate marker positions and REF alleles against the reference with `check_markers_against_fasta.py`.
3. Run `KASP_design_V1.2.py` using the original VCF or a background CSV when available.
4. Keep multiple ranked assays per marker if you want wet-lab choice.
5. For difficult loci with unavoidable shared non-target variants in `A1` and `A2`, consider the flagged shared-variant tolerance mode.
6. Add BLAST screening when specificity is a concern or when you want BLAST-assisted ranking.

---

## Troubleshooting

### Output is empty or unexpectedly small
Check these first:
- FASTA chromosome names really match the marker table after `chr` normalization
- positions are 1-based
- marker rows use only supported SNP or simple indel formats
- the REF allele truly matches the FASTA for indels
- the background VCF or CSV is not so dense that nearly all primer candidates overlap another variant
- if the locus has unavoidable shared non-target variants in both `A1` and `A2`, consider whether you want to enable the flagged tolerance mode
- `--include_fallback` is off by default

### `Output rows written: 0`
This means no marker produced a selected row under the current settings. Inspect the terminal `Outcome summary`, because the script records the dominant failure reason for each marker.

### Many markers fail as `ref_alt_do_not_match_fasta`
Possible causes: wrong reference genome version, positions are not 1-based, chromosome names do not match, or the marker CSV came from a different coordinate system. For indels, also check that your deletion range or anchored REF allele really matches the FASTA span.

### Many markers fail as `equal_length_multibase_not_supported`
This script currently does not support substitutions like `AT -> GC`. Split them into SNPs if biologically appropriate, or extend the code.

### Many markers fail as `*_overlaps_variant`
This usually means the surrounding region is highly polymorphic. You may need to choose a different marker set, reduce background density, or increase the number of candidate regions scanned.

For difficult loci where the **same unavoidable non-target position(s)** overlap both `A1` and `A2`, you can try:
- `--allow_shared_allele_primer_variants`
- `--max_shared_allele_variant_positions 1`
- `--min_shared_allele_variant_distance_3prime 5`

Returned rows that use this exception are clearly flagged in the output.

### Human-readable deletion rows are being skipped
Check:
- the position really uses `start-end`
- the `ref` length matches the range length
- the `alt` is empty, `-`, `.`, `DEL`, or `DELETION`

### The script says raw VCF is not accepted as `--marker_csv`
That is expected in the current version. Keep `--marker_csv` as a CSV/TSV marker table and use the VCF only for:
- conversion with `vcf_to_kasp_csv.py`
- background masking with `--background_vcf`

For indels from a VCF, copy or convert the desired rows into anchored marker-table form such as:

```csv
Chr,position,ref,alt
ChrA01,30423823,TAA,TCTTAAA
```


### BLAST errors
Make sure `blastn` is installed and visible in `PATH`, `makeblastdb` is installed if auto-building is enabled, the FASTA is readable, and the BLAST database path points to the database prefix rather than only a sidecar file.

---

## Design philosophy

This repository is intentionally conservative.

It prefers:
- rejecting primers that overlap known nearby variants by default
- allowing only a narrow, explicitly flagged exception for shared non-target positions in both `A1` and `A2` when you enable it
- retaining multiple ranked designs per marker
- using Primer3 for thermodynamic checks
- optionally adding BLAST-based off-target visibility

rather than simply returning the first sequence pair that appears to work. That makes the output more selective, but generally better suited for downstream assay testing.

---

## Appendix: full output columns in `--output_mode full`

The list below annotates the major columns written when the script is run with `--output_mode full`.

### General notes
- `A1` is the REF-linked assay primer in the current code
- `A2` is the ALT-linked assay primer in the current code
- all genomic coordinates are 1-based
- all `*_tm` values are Primer3-derived temperatures in degrees C
- all `*_gc` values are GC percentages
- `Common_Primer` and `Common_Core` are the same sequence in the current implementation

### Marker identity and normalization

| Column | Meaning |
|---|---|
| `Chr` | Chromosome or contig name from the normalized marker table. |
| `Position` | Internal 1-based coordinate used by the design engine. For deletion ranges, this is the start of the deleted block. |
| `Input_position` | Original input token from the marker CSV, such as `13969`, `1300`, or `13969-13971`. |
| `Ref` | REF allele from the marker CSV after normalization to uppercase. |
| `Alt` | ALT allele from the marker CSV after normalization to uppercase. |
| `FASTA_base` | FASTA base or FASTA slice reported by the code at the target marker region. |
| `FASTA_allele_annotation` | Consistency label reported by the code. |
| `TASSEL_major_allele_as_REF` | `Yes`, `No`, or `Unknown`. |

### Primer sequences and orientation

| Column | Meaning |
|---|---|
| `Orientation` | Assay layout direction. SNP and indel orientations use different labels. |
| `A1_Primer` | Full sequence to order for the REF-linked primer: FAM tail + `A1_Core`. |
| `A2_Primer` | Full sequence to order for the ALT-linked primer: VIC tail + `A2_Core`. |
| `Common_Primer` | Full common primer sequence to order. |
| `A1_Core` | Genomic binding region of A1 without the tail. |
| `A2_Core` | Genomic binding region of A2 without the tail. |
| `Common_Core` | Genomic binding region of the common primer. |

### Primer genomic spans and basic assay metrics

| Column | Meaning |
|---|---|
| `A1_span_start` | 1-based chromosome start coordinate where `A1_Core` binds. |
| `A1_span_end` | 1-based chromosome end coordinate where `A1_Core` binds. |
| `A2_span_start` | 1-based chromosome start coordinate where `A2_Core` binds. |
| `A2_span_end` | 1-based chromosome end coordinate where `A2_Core` binds. |
| `Allele_specific_span_start` | Minimum of `A1_span_start` and `A2_span_start`. |
| `Allele_specific_span_end` | Maximum of `A1_span_end` and `A2_span_end`. |
| `Common_primer_span_start` | 1-based chromosome start coordinate where the common primer binds. |
| `Common_primer_span_end` | 1-based chromosome end coordinate where the common primer binds. |
| `A1_core_tm` | Melting temperature of `A1_Core`. |
| `A2_core_tm` | Melting temperature of `A2_Core`. |
| `Common_core_tm` | Melting temperature of `Common_Core`. |
| `A1_core_gc` | GC percentage of `A1_Core`. |
| `A2_core_gc` | GC percentage of `A2_Core`. |
| `Common_core_gc` | GC percentage of `Common_Core`. |
| `Product_size` | Expected PCR amplicon size in base pairs using the script's current internal calculation. |

### Shared-variant tolerance columns

These columns are populated for all rows and are especially useful when `--allow_shared_allele_primer_variants` is enabled.

| Column | Meaning |
|---|---|
| `Shared_variant_tolerance_used` | `Yes` if the optional shared-variant exception was used for this row, otherwise `No`. |
| `Shared_variant_positions` | Semicolon-separated list of tolerated shared non-target variant positions overlapping both `A1` and `A2`. |
| `Shared_variant_count` | Number of tolerated shared non-target variant positions. |
| `Shared_variant_warning` | Human-readable explanation of the tolerated shared-overlap condition. |
| `A1_shared_variant_min_3prime_distance` | Minimum distance in nucleotides from the 3' end of `A1` to any tolerated shared non-target variant position. |
| `A2_shared_variant_min_3prime_distance` | Minimum distance in nucleotides from the 3' end of `A2` to any tolerated shared non-target variant position. |

### Thermodynamic and structure-check columns

| Column | Meaning |
|---|---|
| `passed` | `True` if the assay met strict thermodynamic and size parameters. |
| `Exclusion_reason` | Reason the row failed local PASS criteria when present. |
| `worst_struct_tm` | Highest structure Tm observed across checked hairpins, homodimers, and heterodimers. |
| `max_tm_delta` | Largest absolute Tm difference among `A1_Core`, `A2_Core`, and `Common_Core`. |
| `A1_hairpin_tm` | Primer3 predicted hairpin Tm for `A1_Primer`. |
| `A2_hairpin_tm` | Primer3 predicted hairpin Tm for `A2_Primer`. |
| `Common_hairpin_tm` | Primer3 predicted hairpin Tm for `Common_Primer`. |
| `A1_homodimer_tm` | Primer3 predicted homodimer Tm for `A1_Primer`. |
| `A2_homodimer_tm` | Primer3 predicted homodimer Tm for `A2_Primer`. |
| `Common_homodimer_tm` | Primer3 predicted homodimer Tm for `Common_Primer`. |
| `A1_Common_heterodimer_tm` | Primer3 predicted heterodimer Tm between `A1_Primer` and `Common_Primer`. |
| `A2_Common_heterodimer_tm` | Primer3 predicted heterodimer Tm between `A2_Primer` and `Common_Primer`. |
| `A1_A2_heterodimer_tm` | Primer3 predicted heterodimer Tm between `A1_Primer` and `A2_Primer`. |

### Local ranking and final selection columns

| Column | Meaning |
|---|---|
| `Local_score` | Weighted penalty score computed from structure Tm, Tm balance, product size, and GC. Lower is better. |
| `Local_score_tuple` | Pipe-delimited raw penalty metrics used for internal sorting. |
| `Ranking_mode` | Logic used to rank rows, such as `local_only` or `local_plus_blast`. |
| `Rank` | 1-based rank among the rows written for that marker. |
| `Combined_rank_key` | Final sort key written as text. |
| `Design_status` | `PASS` or `FALLBACK`. |

### BLAST-only off-target columns

These columns are added only when BLAST screening is enabled with `--blast_db` or `--blast_use_genome`.

| Column | Meaning |
|---|---|
| `A1_off_target_3prime_hits` | Number of off-target matches covering the 3' end of `A1_Core`. |
| `A2_off_target_3prime_hits` | Number of off-target matches covering the 3' end of `A2_Core`. |
| `Common_off_target_3prime_hits` | Number of off-target matches covering the 3' end of `Common_Core`. |
| `A1_exact_full_length_off_target_hits` | Number of exact full-length off-target matches for `A1_Core`. |
| `A2_exact_full_length_off_target_hits` | Number of exact full-length off-target matches for `A2_Core`. |
| `Common_exact_full_length_off_target_hits` | Number of exact full-length off-target matches for `Common_Core`. |
| `A1_best_offtarget_3prime_bitscore` | Highest BLAST bitscore among 3'-covering off-targets for `A1_Core`. |
| `A2_best_offtarget_3prime_bitscore` | Highest BLAST bitscore among 3'-covering off-targets for `A2_Core`. |
| `Common_best_offtarget_3prime_bitscore` | Highest BLAST bitscore among 3'-covering off-targets for `Common_Core`. |
| `Total_off_target_3prime_hits` | Sum of all 3'-covering hits across all three primers. |
| `Total_exact_full_length_off_target_hits` | Sum of all exact full-length off-target matches across all three primers. |
| `Worst_offtarget_3prime_bitscore` | Maximum bitscore among all 3'-covering off-targets across all three primers. |
