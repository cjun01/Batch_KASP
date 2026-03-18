
# KASP Marker Design Tool

## Overview
The KASP Marker Design Tool is a Python-based utility for designing KASP (Kompetitive Allele Specific PCR) primers from SNP markers. It includes:

- `vcf_to_kasp_csv.py` to convert a VCF into the marker CSV format required by the designer.
- `KASP_design.py` to design FAM/VIC allele-specific forward primers plus a common reverse primer.

## Features
- Convert VCF or VCF.gz files into KASP marker input tables.
- Automated primer design for KASP markers.
- Primer-region background variant masking to avoid nearby known variants inside primer sites.
- Simple Tm, GC, repeat, hairpin, and primer-dimer checks.

## Requirements
- Python 3.10.13 or above
- Biopython
- Pandas

## Installation
1. Ensure Python 3 is installed.
2. Install the required packages:

```bash
pip install -r requirements.txt
```

## Workflow
The typical workflow is:

1. Convert a VCF into a marker CSV with `vcf_to_kasp_csv.py`.
2. Select a target marker list for assay design.
3. Run `KASP_design.py` on that target marker CSV and a reference genome FASTA.
4. Optionally provide the original VCF as background variation so primer sites can avoid nearby known variants.

## Quick start
Run these commands from the `Batch_KASP` folder.

### 1. Create a marker CSV from a VCF
```bash
python vcf_to_kasp_csv.py \
  -i input.vcf.gz \
  -o all_markers.csv
```

### 2. Prepare a target list
`KASP_design.py` works best when `-m/--marker_csv` is a curated target list of SNPs you actually want assays for.

Do not use a whole-genome SNP list as the target list unless that is really your intention. If you give the script hundreds of thousands of SNPs as targets, many or all of them may be rejected because nearby variants fall inside candidate primer regions.

Your target CSV must contain:

- `Chr`
- `position`
- `ref`
- `alt`

### 2b. Optional: check the target list against the FASTA first
If you want to verify that the target coordinates and `ref` alleles match the reference genome before primer design, run:

```bash
python check_markers_against_fasta.py \
  -g genome.fa \
  -m target_markers.csv \
  -o marker_fasta_check.csv
```

This prints a summary and writes a detailed report with per-marker status such as `match`, `ref_mismatch`, `missing_chromosome_in_fasta`, or `position_out_of_range`.

### 3. Run primer design
Recommended command:

```bash
python KASP_design.py \
  -g genome.fa \
  -m target_markers.csv \
  -o primer_designs.csv \
  --background_vcf input.vcf.gz
```

This uses:

- `target_markers.csv` as the SNPs you want assays for
- `input.vcf.gz` as the background variation to avoid inside primer binding regions

### 4. Check the summary at the end
The script prints a summary like:

```text
Designed X primer sets from Y targets.
Skipped Z targets with no acceptable primer pair.
Reference mismatches against FASTA: N
Outcome summary:
  designed: ...
  forward_primer_overlaps_variant: ...
  reverse_primer_overlaps_variant: ...
  missing_chromosome_in_fasta: ...
```

If you get `Designed 0 primer sets`, the summary is the first thing to inspect.

## Step 1: Convert a VCF to KASP marker CSV
`KASP_design.py` expects a CSV with these columns:

- `Chr`
- `position`
- `ref`
- `alt`

The helper script `vcf_to_kasp_csv.py` creates that file directly from a VCF.

### Basic conversion
Use this when you want to convert the whole VCF into KASP input markers.
It keeps all simple biallelic SNPs in the file and skips indels or multiallelic sites.

```bash
python vcf_to_kasp_csv.py \
  -i input.vcf.gz \
  -o markers_for_kasp.csv
```

### Sample-specific conversion
Use this when the VCF contains many samples but you want markers for one sample only.
This command checks the genotype of `SampleName` and keeps only variants that pass the sample filter.

```bash
python vcf_to_kasp_csv.py \
  -i input.vcf.gz \
  -o sample_markers.csv \
  -s SampleName
```

By default, when `--sample` is provided, the script keeps only homozygous alternate (`1/1`) SNPs for that sample.

In short, the difference between the two commands is:

- Without `-s/--sample`: keep all simple SNPs from the VCF.
- With `-s/--sample`: keep only SNPs relevant to the named sample.

### Optional filters
- `--genotype alt`: Keep only `1/1` calls for the selected sample. This is the default when `--sample` is used.
- `--genotype nonref`: Keep `0/1`, `1/0`, and `1/1` calls.
- `--genotype any`: Keep any non-missing genotype call.
- `--exclude-contig-substring unitig`: Skip contigs whose names contain a substring such as `unitig`.

Example:

```bash
python vcf_to_kasp_csv.py \
  -i input.vcf.gz \
  -o sample_markers.csv \
  -s SampleName \
  --genotype nonref \
  --exclude-contig-substring unitig
```

### Output format from `vcf_to_kasp_csv.py`

```csv
Chr,position,ref,alt
Lcu.2RBY.Chr1,95213,A,G
Lcu.2RBY.Chr2,128004,C,T
```

## Step 2: Run KASP primer design
Run `KASP_design.py` with a reference genome FASTA and the CSV created in step 1.

```bash
python KASP_design.py \
  -g genome.fa \
  -m markers_for_kasp.csv \
  -o primer_designs.csv
```

### Recommended primer-safe run
This is the recommended mode when you have the original VCF available. The marker CSV gives the target SNPs, and the background VCF is used to reject primer candidates that overlap other known variants in the primer binding regions.

```bash
python KASP_design.py \
  -g genome.fa \
  -m sample_markers.csv \
  -o primer_designs.csv \
  --background_vcf input.vcf.gz
```

### Command line arguments for `KASP_design.py`
- `-g` or `--genome_path`: Path to the reference genome FASTA file.
- `-m` or `--marker_csv`: Path to the target marker CSV with columns `Chr,position,ref,alt`.
- `-o` or `--output_csv`: Path to save the primer design results.
- `--background_vcf`: Optional VCF or VCF.gz of background variants to avoid inside primer binding regions.
- `--background_csv`: Optional CSV of background variants to avoid inside primer binding regions.

### How nearby SNP handling works now
`KASP_design.py` no longer relies only on a fixed "previous SNP / next SNP" distance rule. Instead, it:

- uses the target marker CSV as a minimum background set
- optionally adds all variants from `--background_vcf` or `--background_csv`
- rejects candidate primer pairs if another known variant overlaps either primer binding region

This is more biologically meaningful than filtering only by marker-to-marker distance, because a nearby SNP matters mainly when it lands inside a primer.

### Important input rules
- The chromosome names in the marker CSV should refer to the same chromosomes as the FASTA. The script tolerates an optional `chr` prefix, so `1LG6` and `chr1LG6` are treated as the same chromosome.
- The `position` column must be 1-based genomic coordinates.
- The `ref` allele should match the reference FASTA at that position.
- The `-m/--marker_csv` file should usually be a target list, not every SNP from the VCF.
- Use `--background_vcf` for `.vcf` or `.vcf.gz` files and `--background_csv` for `.csv` files.

If chromosome names do not match, all targets may be skipped with `missing_chromosome_in_fasta`.

## Output file format
The output CSV from `KASP_design.py` contains:

- `Chr`: Chromosome.
- `Position`: SNP position.
- `Ref`: Reference allele.
- `Alt`: Alternate allele.
- `A1_Primer`: Forward primer labeled with FAM.
- `A2_Primer`: Forward primer labeled with VIC.
- `GC content for forward primer`: GC content of the forward primer.
- `Repeats within forward primer`: Repeat check result for the forward primer.
- `Reverse_Primer`: Common reverse primer.
- `GC content for reverse primer`: GC content of the reverse primer.
- `Repeats within reverse primer`: Repeat check result for the reverse primer.
- `Tm`: Primer melting temperature.
- `Product size`: Expected amplicon size.
- `Hairpin`: Hairpin check result.
- `Primer dimer`: Primer dimer check result.
- `SNP alignment Check`: Whether the reference base in the FASTA matches the provided `ref` allele.

## End-to-end example

```bash
python vcf_to_kasp_csv.py \
  -i input.vcf.gz \
  -o sample_markers.csv \
  -s SampleName \
  --exclude-contig-substring unitig

python KASP_design.py \
  -g genome.fa \
  -m sample_markers.csv \
  -o primer_designs.csv \
  --background_vcf input.vcf.gz
```

## Troubleshooting
### `Designed 0 primer sets from N targets`
Common causes:

- `-m/--marker_csv` contains all SNPs from the VCF instead of a smaller target list
- chromosome names in the marker CSV do not match the FASTA
- the `ref` allele in the marker CSV does not match the FASTA
- nearby variants from the target list or background VCF overlap candidate primer regions

Check the printed `Outcome summary` to see the main failure reason.
