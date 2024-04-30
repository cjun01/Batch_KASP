
# KASP Marker Design Tool

## Overview
The KASP Marker Design Tool is a Python-based utility designed to assist in the design of KASP (Kompetitive Allele Specific PCR) markers. It automates the process of designing forward (FAM and VIC labeled) and reverse primers for a given set of genetic markers, taking into consideration factors like melting temperature (Tm), potential hairpin structures, and possible primer dimer.

## Features
- Automated primer design for KASP markers.
- Customizable Tm calculation and hairpin structure and primer dimer checking.
- Command-line interface for easy use and integration into bioinformatics pipelines.

## Requirements
- Python 3.10.13 or above
- Biopython
- Pandas

## Installation
1. Ensure that Python 3.x is installed on your system.
2. Install required Python packages.

## Usage
The tool can be executed from the command line by providing paths to the input genome FASTA file, the CSV file containing marker information, and the desired path for the output CSV file containing the primer designs.

### Command Line Arguments
- `-g` or `--genome_path`: Path to the genome FASTA file.
- `-m` or `--marker_csv`: Path to the CSV file containing marker information. This file should have columns for Chromosome (`Chr`), Position (`Position`), and Alternate allele (`Alt`).
- `-o` or `--output_csv`: Path to save the output CSV file with primer designs.

### Running the Tool
## Input File Format
The marker information CSV file should contain three columns:
- `Chr`: Chromosome name or number.
- `Position`: Position of the SNP.
- `Alt`: Alternate allele at SNP position.

Example of marker CSV: SNP_position.csv

## Output File Format (primer_designs.csv)
The output CSV file will contain the following columns:
- `Chr`: Chromosome.
- `Position`: Position of the SNP.
- `Alt`: Alternate allele.
- `A1_Primer`: Forward primer labeled with FAM.
- `A2_Primer`: Forward primer labeled with VIC, specific to the alternate allele.
- `Reverse_Primer`: Common reverse primer.
- `Tm`: Melting temperature primers in PCR.
- `Product size`: Expected size of the amplified product.
- `hairpin`: Check for potential hairpin structure.
- `Primer dimer`: Check for the possible formation of primer dimer.
