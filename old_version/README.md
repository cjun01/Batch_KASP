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

1.  **Loads the reference genome**
    - accepts FASTA or FASTA.gz
    - normalizes chromosome names so `chr1` and `1` are treated as the same chromosome

2.  **Reads and cleans the marker table**
    - requires `Chr`, `position`, `ref`, and `alt`
    - treats column names case-insensitively
    - keeps only simple biallelic A/C/G/T SNPs
    - removes exact duplicate marker rows after normalization

3.  **Checks allele consistency against the FASTA**
    - `REF matches FASTA`
    - `ALT matches FASTA; FASTA carries ALT allele`
    - `Neither REF nor ALT matches FASTA`
    
    Markers in the third category are rejected.

4.  **Builds a strict variant mask for primer placement**
    - all target SNPs are masked automatically
    - optional background variants from a VCF or CSV are also masked
    - the current target SNP itself is still allowed inside its own allele-specific primer region
    - any other overlapping known variant causes that primer candidate to be rejected

5.  **Searches both supported KASP assay layouts**
    - **Orientation 1**: allele-specific forward primers + common reverse primer
    - **Orientation 2**: allele-specific reverse primers + common forward primer

6.  **Adds the standard reporter tails**
    - `A1_Primer = FAM tail + REF-specific core`
    - `A2_Primer = VIC tail + ALT-specific core`
    
    Default tails:
    - FAM tail: `GAAGGTGACCAAGTTCATGCT`
    - VIC tail: `GAAGGTCGGAGTCAACGGATT`

7.  **Evaluates thermodynamic quality with Primer3**
    - hairpins
    - homodimers
    - heterodimers
    - Tm balance among primer cores

8.  **Ranks local candidates**
    - candidates are scored and sorted using a strict local ranking rule
    - multiple ranked candidates can be retained per SNP

9.  **Optionally performs BLAST off-target screening**
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

### 2) Marker CSV

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

Must contain at least `Chr` and `position`. If a `ref` column is present, the masked span uses its length. Otherwise, a single base is masked.

**Important masking rule**

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
python vcf_to_kasp_csv.py   -i input.vcf.gz   -o target_markers.csv
```

Example: sample-specific conversion
Keep only markers for a specific sample:

```bash
python vcf_to_kasp_csv.py   -i input.vcf.gz   -o target_markers.csv   -s SampleName
```

Keep only non-reference genotypes for that sample:

```bash
python vcf_to_kasp_csv.py   -i input.vcf.gz   -o target_markers.csv   -s SampleName   --genotype nonref
```

---

## Step 2: Check markers against the FASTA

Before design, verify that the marker coordinates and `REF` alleles are consistent with the reference.

```bash
python check_markers_against_fasta.py   -g genome.fa   -m target_markers.csv   -o marker_fasta_check.csv
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
python KASP_design.py   -g genome.fa   -m target_markers.csv   -o primer_designs.csv
```

By default, the current script now writes a **concise output table**: it keeps the strict internal screening and ranking logic, but only writes a smaller set of decision-useful columns. Use `--output_mode full` if you want the full diagnostic table.

If a marker does not receive any written assay row, the script also writes a companion failed-marker report by default using the output stem plus `_failed.csv`. For example, `primer_designs.csv` is paired with `primer_designs_failed.csv`.

### Example: include a background VCF (Recommended)

This is usually the recommended mode when your markers came from a larger variant set.

```bash
python KASP_design.py   -g genome.fa   -m target_markers.csv   -o primer_designs.csv   --background_vcf input.vcf.gz
```

### Example: stricter local filters

```bash
python KASP_design.py   -g genome.fa   -m target_markers.csv   -o primer_designs.csv   --background_vcf input.vcf.gz   --product_min 60   --product_max 120   --max_thermo_tm 40
```

This configuration narrows the acceptable amplicon size and lowers tolerance for strong unwanted structures.

---

## Step 4: Add BLAST off-target screening

BLAST screening is optional and is mainly useful when the genome is large or repetitive, duplicated loci are a concern, or you want BLAST-assisted ranking among already acceptable local candidates.

### Option A: use an existing BLAST database

```bash
python KASP_design.py   -g genome.fa   -m target_markers.csv   -o primer_designs.csv   --background_vcf input.vcf.gz   --blast_db /path/to/blast_db_prefix
```

### Option B: use the genome FASTA directly and auto-build the BLAST database

```bash
python KASP_design.py   -g genome.fa   -m target_markers.csv   -o primer_designs.csv   --background_vcf input.vcf.gz   --blast_use_genome
```

---

## Candidate generation and ranking

The script does **not** stop at the first acceptable assay. Instead, it can keep multiple acceptable candidates per SNP and rank them.

### Local candidate generation

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

The script writes a scalar `Local_score`, but the primary local ordering is driven by this tuple. This allows the pipeline to return several acceptable designs for the same SNP and rank them from better to worse.

---

## `PASS` vs `FALLBACK`

Each candidate is first assessed locally. The boolean field `passed` is `True` only if both conditions are satisfied:
- `product_min <= Product_size <= product_max`
- `worst_struct_tm <= max_thermo_tm`

Selected output rows use `Design_status`:
- `PASS` — selected from candidates that satisfy the local pass criteria
- `FALLBACK` — selected from weaker candidates only when no `PASS` candidates exist

By default, only `PASS` rows are written. To include fallback rows when a marker has no passing candidate:

```bash
python KASP_design.py   -g genome.fa   -m target_markers.csv   -o primer_designs.csv   --include_fallback
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

**Output retention settings**
- `--top_n` — default `5`; maximum number of written rows per marker
- `--include_fallback` — write fallback rows when no `PASS` rows exist for that marker
- `--preblast_top_n` — default `10`; with BLAST enabled, only this many best local candidates per marker are BLAST-screened
- `--output_mode` — `concise` or `full`; default `concise`
- `--failed_output_csv` — optional path for the companion failed-marker report

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

By default, `KASP_design.py` now writes a **concise** output schema. This keeps the design logic just as strict, but reduces the number of reported indicators after the main assay-definition and basic metric fields. Use `--output_mode full` if you want the full diagnostic table described in the Appendix.

The output columns change dynamically depending on whether the assay design was successful, whether you enabled BLAST off-target screening, and whether you chose the "concise" or "full" output mode.

Below is a comprehensive guide to what the default concise columns mean, grouped by category. *(For the full list of extended diagnostic columns, see the Appendix).*

### 1. Genomic Context (The Target SNP)

These columns describe the biological variant you are targeting and how it maps to your reference genome.

- **`Chr`**: The chromosome or contig where the target SNP is located.
- **`Position`**: The 1-based genomic coordinate of the target SNP.
- **`Ref`**: The reference allele (nucleotide) provided in your input marker CSV.
- **`Alt`**: The alternate allele provided in your input marker CSV.
- **`FASTA_base`**: The actual nucleotide base found at this exact position in the provided reference genome (`.fasta`).
- **`FASTA_allele_annotation`**: A diagnostic message indicating whether the genome sequence carries the `Ref` allele, the `Alt` allele, or neither. (Crucial for catching reference assembly or strand issues).
- **`TASSEL_major_allele_as_REF`**: Outputs "Yes", "No", or "Unknown". Bioinformatics pipelines (like TASSEL) sometimes force the most common allele in a population into the "Ref" column, even if the physical reference genome carries the "Alt". This column flags `"Yes"` if your `Alt` allele is the one actually found in the FASTA genome.

### 2. Primer Sequences to Order

KASP requires three primers: two allele-specific forward (or reverse) primers, and one common reverse (or forward) primer.

- **`Orientation`**: The physical direction of the assay on the DNA strand.
    - *`allele_specific_forward_common_reverse`*: The two allele primers bind to the top strand; the common primer binds downstream on the bottom strand.
    - *`allele_specific_reverse_common_forward`*: The two allele primers bind to the bottom strand; the common primer binds upstream on the top strand.
- **`A1_Primer`**: The full sequence to order for Allele 1 (`Ref`). It includes the genomic binding sequence **plus the standard FAM tail** (`GAAGGTGACCAAGTTCATGCT`) at the 5' end.
- **`A2_Primer`**: The full sequence to order for Allele 2 (`Alt`). It includes the genomic binding sequence **plus the standard VIC tail** (`GAAGGTCGGAGTCAACGGATT`) at the 5' end.
- **`Common_Primer`**: The full sequence to order for the common flanking primer. It does not have a tail.
- **`A1_Core` / `A2_Core` / `Common_Core`**: The exact template-binding sequence of the primers *without* the FAM/VIC tails. This is the sequence actually used to calculate thermodynamics and genomic coordinates.

### 3. Genomic Coordinates & Product Size

- **`Allele_specific_span_start` / `Allele_specific_span_end`**: The 1-based start and end coordinates on the chromosome where the allele-specific core primers bind. (The target SNP is located at the 3' end of this span).
- **`Common_primer_span_start` / `Common_primer_span_end`**: The 1-based start and end coordinates on the chromosome where the common primer binds.
- **`Product_size`**: The expected length (in base pairs) of the final PCR amplicon, measured from the 5' end of the allele core to the 5' end of the common core.

### 4. Thermodynamics (Primer3 Metrics)

These columns assess the physical likelihood that the PCR will work effectively without forming primer-dimers.

- **`A1_core_tm` / `A2_core_tm` / `Common_core_tm`**: The predicted Melting Temperature (Tm) of the core sequences in °C.
- **`A1_core_gc` / `A2_core_gc` / `Common_core_gc`**: The GC content percentage of the core sequences. Extreme GC content causes poor amplification.
- **`worst_struct_tm`**: The script tests 9 different potential secondary structures (e.g., A1 forming a hairpin, A1 binding to A2, A2 binding to the Common primer). This column reports the **highest (worst) Tm** among all of them. A high number means the primers will stick to each other instead of your DNA.
- **`max_tm_delta`**: The maximum temperature difference between the A1, A2, and Common primer Tms. (Lower is better; a large gap means one primer will bind poorly at the shared annealing temperature).

### 5. Scoring & Ranking

These columns explain how the script evaluated the assay compared to other possible designs for the exact same SNP.

- **`passed`**: `True` if the assay met all strict thermodynamic, GC, and size parameters. `False` if it violated a parameter but was kept as a fallback.
- **`Exclusion_reason`**: If `passed` is `False`, this explains why (e.g., `product_size_below_min`).
- **`Local_score`**: A custom penalty score calculated by the script. It adds penalty points for deviating from the ideal Tm/GC/Product size, and heavily penalizes secondary structures and Tm gaps. **Lower scores are better.** (0.0 is perfect).
- **`Ranking_mode`**: Explains the logic used to rank the primers (e.g., `local_only` based purely on thermodynamics, or `local_plus_blast` incorporating genome specificity).
- **`Rank`**: The final rank of this primer set for the specific SNP (1 is the top recommendation).
- **`Design_status`**:
    - `PASS`: The assay met all your strict limits.
    - `FALLBACK`: Failed strict thresholds, but was provided as the "best available option" because you ran the script with `--include_fallback`.

### 6. BLAST Off-Target Metrics (Only if BLAST is enabled)

If you ran the script with `--blast_use_genome` or `--blast_db`, these columns warn you if the primers will amplify the wrong parts of the genome.

- **`Total_exact_full_length_off_target_hits`**: The number of times the primers perfectly match another sequence in the genome from end to end. (Ideally 0).
- **`Total_off_target_3prime_hits`**: The total number of off-target matches that cover the **3' end** of your primers. Because Taq polymerase extends from the 3' end, a 3' match is highly dangerous and likely to cause accidental background PCR amplification.
- **`Worst_offtarget_3prime_bitscore`**: The highest BLAST bitscore among all the 3' off-target hits. A high score means a very strong, stable mis-priming event is possible.

---

## Companion failed-marker report

Markers that do not receive any written assay rows (even as a fallback) are written to a separate CSV. By default, the path is derived from the main output path:

- main output: `primer_designs.csv`
- failed-marker report: `primer_designs_failed.csv`

The failed-marker report is highly useful for diagnosing exactly why a genomic region is hostile to primer design. It includes:

- **`Chr`, `Position`, `Ref`, `Alt`**: The original marker details.
- **`FASTA_base`, `FASTA_allele_annotation`, `TASSEL_major_allele_as_REF`**: The FASTA alignment check results.
- **`Failure_reason`**: The primary reason the script gave up (e.g., `ref_alt_do_not_match_fasta`, `position_out_of_range`, or `no_candidate_found`).
- **`Failure_reason_counts`**: A highly detailed tally of why candidate primers were rejected during the internal scanning loop. For example: `forward_allele_specific_tm=12; forward_allele_specific_overlaps_variant=5`. This tells you *exactly* why the region failed (e.g., too AT-rich, or too many background SNPs nearby blocking primer placement).
- **`PASS_candidates_found` / `FALLBACK_candidates_found`**: Indicates how many candidates were successfully generated. (Usually `0` if the marker ended up in this file, but helps diagnose if strict filtering removed them later).

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

*(Note: `Reference mismatches against FASTA` is counted across all processed markers, not only written output rows).*

---

## Practical examples

### Example 1: simplest possible run

```bash
python KASP_design.py -g genome.fa -m markers.csv -o kasp_out.csv
```

Use this when you already trust the marker file and do not need background masking or BLAST.

### Example 2: recommended design run with variant masking

```bash
python KASP_design.py -g genome.fa -m markers.csv -o kasp_out.csv --background_vcf all_variants.vcf.gz
```

Use this when your markers were selected from a larger variant call set and you want to avoid adjacent known variants inside primer-binding regions.

### Example 3: stricter amplicon and structure filtering

```bash
python KASP_design.py -g genome.fa -m markers.csv -o kasp_out.csv   --background_vcf all_variants.vcf.gz   --product_min 70   --product_max 110   --max_thermo_tm 38
```

### Example 4: keep more candidate rows per marker

```bash
python KASP_design.py -g genome.fa -m markers.csv -o kasp_out.csv   --background_vcf all_variants.vcf.gz   --top_n 10
```

### Example 5: include fallback rows

```bash
python KASP_design.py -g genome.fa -m markers.csv -o kasp_out.csv   --background_vcf all_variants.vcf.gz   --include_fallback
```

### Example 6: BLAST-assisted ranking with auto database creation

```bash
python KASP_design.py -g genome.fa -m markers.csv -o kasp_out.csv   --background_vcf all_variants.vcf.gz   --blast_use_genome   --top_n 5   --product_min 60   --product_max 120   --max_thermo_tm 40   --preblast_top_n 10   --blast_num_threads 4
```

### A Real Case

```bash
python /home/cflzxc/software/Batch_KASP/KASP_design.py   -g /output/genomic/plant/Pisum/sativum/Genome/Pisum_sativum_v1a.fa   -m markers_for_kasp.csv   -o primer_designs.csv   --background_csv ALL_SNPs.csv   --blast_use_genome   --top_n 5   --product_min 60   --product_max 120   --max_thermo_tm 40   --preblast_top_n 10   --blast_num_threads 4
```

**One-liner (with fallbacks included):**

```bash
python /home/cflzxc/software/Batch_KASP/KASP_design.py -g /output/genomic/plant/Pisum/sativum/Genome/Pisum_sativum_v1a.fa -m markers_for_kasp.csv -o primer_designs.csv --background_csv ALL_SNPs.csv --blast_use_genome --top_n 5 --product_min 60 --product_max 120 --max_thermo_tm 40 --preblast_top_n 10 --blast_num_threads 4 --include_fallback
```

---

## Recommended workflow

1. Convert your VCF to a clean marker CSV with `vcf_to_kasp_csv.py`.
2. Validate marker positions and REF alleles against the reference with `check_markers_against_fasta.py`.
3. Run `KASP_design.py` using the original VCF as `--background_vcf` when available.
4. Keep multiple ranked assays per marker if you want choice for wet-lab testing.
5. Add BLAST screening when specificity is a concern or when you want BLAST-assisted ranking.

---

## Troubleshooting

### Output is empty or unexpectedly small
Check these first:
- FASTA chromosome names really match the marker table after `chr` normalization
- marker positions are 1-based
- the marker file contains only simple SNPs, not indels or multiallelic rows
- at least one of `REF` or `ALT` matches the FASTA
- the background VCF or CSV is not so dense that nearly all primer candidates overlap another variant
- `--include_fallback` is off by default

### `Output rows written: 0`
This means no marker produced a selected row under the current settings. Inspect the terminal `Outcome summary`, because the script records the dominant failure reason for each marker.

### Many markers fail as `ref_alt_do_not_match_fasta`
Possible causes: wrong reference genome version, positions are not 1-based, chromosome names do not match, or marker CSV came from a different coordinate system. Run `check_markers_against_fasta.py` to confirm.

### Many markers fail as `*_overlaps_variant`
This usually means the surrounding region is highly polymorphic. You may need to choose a different SNP set, reduce background density, or increase the number of candidate regions scanned.

### BLAST errors
Make sure `blastn` is installed and visible in `PATH`, `makeblastdb` is installed if auto-building is enabled, the FASTA is readable, and the BLAST database path points to the database prefix (not only a single sidecar file).

---

## Design philosophy

This repository is intentionally conservative.

It prefers:
- rejecting primers that overlap known nearby variants
- retaining multiple ranked designs per SNP
- using Primer3 for thermodynamic checks
- optionally adding BLAST-based off-target visibility

rather than simply returning the first sequence pair that appears to work. That makes the output more selective, but generally better suited for downstream assay testing.

---

## Appendix: full output columns from `KASP_design.py`

The list below annotates **all** the columns written when `KASP_design.py` is run with `--output_mode full`. If you run in the default `concise` mode, only a curated subset of these will appear.

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
| `FASTA_base` | Exact nucleotide found in the reference FASTA at `Chr:Position`. |
| `FASTA_allele_annotation` | Consistency label (`REF matches FASTA`, `ALT matches FASTA; FASTA carries ALT allele`, or `Neither`). |
| `TASSEL_major_allele_as_REF` | `Yes` if the FASTA matched your `Alt` allele (common in TASSEL exports), `No` if it matched `Ref`, or `Unknown`. |
| `Orientation` | Assay layout direction (`allele_specific_forward_common_reverse` or `allele_specific_reverse_common_forward`). |
| `A1_Primer` | Full sequence to order for Allele 1 (FAM tail + `A1_Core`). |
| `A2_Primer` | Full sequence to order for Allele 2 (VIC tail + `A2_Core`). |
| `Common_Primer` | Full common primer sequence to order (no tail). |
| `A1_Core` | Genomic binding region of A1 (no tail). |
| `A2_Core` | Genomic binding region of A2 (no tail). |
| `Common_Core` | Genomic binding region of the common primer. |

### Primer genomic spans and basic assay metrics

| Column | Meaning |
|---|---|
| `Allele_specific_span_start` | 1-based chromosome start coordinate where the allele-specific primer binds. |
| `Allele_specific_span_end` | 1-based chromosome end coordinate where the allele-specific primer binds. |
| `Common_primer_span_start` | 1-based chromosome start coordinate where the common primer binds. |
| `Common_primer_span_end` | 1-based chromosome end coordinate where the common primer binds. |
| `A1_core_tm` | Melting temperature of `A1_Core`. |
| `A2_core_tm` | Melting temperature of `A2_Core`. |
| `Common_core_tm` | Melting temperature of `Common_Core`. |
| `A1_core_gc` | GC percentage of `A1_Core`. |
| `A2_core_gc` | GC percentage of `A2_Core`. |
| `Common_core_gc` | GC percentage of `Common_Core`. |
| `Product_size` | Expected PCR amplicon size in base pairs. |

### Thermodynamic and structure-check columns

These columns assess the physical likelihood that the PCR will work effectively without forming primer-dimers.

| Column | Meaning |
|---|---|
| `passed` | `True` if the assay met strict thermodynamic/size parameters. `False` if kept as a fallback. |
| `worst_struct_tm` | Highest (worst) Tm observed across all checked hairpins, homodimers, and heterodimers for this primer set. |
| `max_tm_delta` | Largest absolute Tm difference among `A1_Core`, `A2_Core`, and `Common_Core`. |
| `A1_hairpin_tm` | Primer3 predicted Hairpin Tm for `A1_Primer`. |
| `A2_hairpin_tm` | Primer3 predicted Hairpin Tm for `A2_Primer`. |
| `Common_hairpin_tm` | Primer3 predicted Hairpin Tm for `Common_Primer`. |
| `A1_homodimer_tm` | Primer3 predicted Homodimer Tm for `A1_Primer`. |
| `A2_homodimer_tm` | Primer3 predicted Homodimer Tm for `A2_Primer`. |
| `Common_homodimer_tm` | Primer3 predicted Homodimer Tm for `Common_Primer`. |
| `A1_Common_heterodimer_tm` | Primer3 predicted Heterodimer Tm between `A1_Primer` and `Common_Primer`. |
| `A2_Common_heterodimer_tm` | Primer3 predicted Heterodimer Tm between `A2_Primer` and `Common_Primer`. |
| `A1_A2_heterodimer_tm` | Primer3 predicted Heterodimer Tm between `A1_Primer` and `A2_Primer`. |

### Local ranking and final selection columns

| Column | Meaning |
|---|---|
| `Local_score` | Weighted penalty score computed from structure Tm, Tm balance, product size, and GC. Lower is better. |
| `Local_score_tuple` | Pipe-delimited raw penalty metrics used for internal sorting (`worst_struct_tm|max_tm_delta|product_size_distance|gc_penalty|tm_center_penalty`). |
| `Ranking_mode` | Logic used to rank rows (e.g., `local_only`, `local_plus_blast`). |
| `Rank` | 1-based rank among the rows written for that SNP (1 is best). |
| `Combined_rank_key` | The final sort key written as text (combines BLAST hits and local scores). |
| `Design_status` | `PASS` (met all criteria) or `FALLBACK` (failed but kept via `--include_fallback`). |

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
| `Worst_offtarget_3prime_bitscore` | Maximum bitscore among all 3'-covering off-targets across all three primers. A high score indicates a stable mis-priming event. |
