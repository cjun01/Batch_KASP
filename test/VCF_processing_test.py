import pandas as pd
import os

def read_vcf(path):
    # Read the lines that do not start with '##' (metadata)
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    # Use pandas to read the lines that are not metadata (header + data)
    return pd.read_csv(pd.io.common.StringIO(''.join(lines)), sep='\t', dtype=str)

def transform_genotype(genotype):
    if genotype == '0/0':
        return 'Ref'
    elif genotype == '1/1':
        return 'Alt'
    elif genotype == '0/1' or genotype == '1/0':
        return 'H'
    elif genotype == './.':
        return 'Missing'  # Added this line to explicitly handle missing genotype calls
    else:
        return genotype  # Handle other cases or non-standard codes

# Usage
vcf_path = r'C:\genome\LR-66_marker_design\LR-66-vcf.txt'
df = read_vcf(vcf_path)
sample_name = 'LR-66-590'  # Replace with your specific sample ID
sample_index = df.columns.get_loc(sample_name)

# Extract genotype information and transform it
df[sample_name] = df.apply(lambda row: transform_genotype(row[sample_index].split(':')[0]), axis=1)
# Filter out rows where the genotype is 'Missing' or 'H'
df = df[~df[sample_name].isin(['Missing', 'H','Ref'])]
# Further filter out rows where the chromosome column contains 'unitig'
df = df[~df['#CHROM'].str.contains('unitig')]
# Remove rows where the reference allele has more than one nucleotide
df = df[df['REF'].str.len() == 1]
# Select relevant columns (example: CHROM, POS, ID, REF, ALT, and the transformed genotype)
df.rename(columns={'#CHROM': 'Chr', 'POS': 'position'}, inplace=True)
output_df = df[['Chr', 'position', 'REF', 'ALT', sample_name]]

# Save to CSV
output_csv_path = r'C:\genome\LR-66_marker_design\LR_66_590_SNP_position.csv'
output_df.to_csv(output_csv_path, index=False)
