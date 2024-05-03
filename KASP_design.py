import argparse
from Bio import SeqIO
import pandas as pd
from tqdm import tqdm  # Import tqdm for the progress bar
FAM = 'GAAGGTGACCAAGTTCATGCT'
VIC = 'GAAGGTCGGAGTCAACGGATT'
class MarkerDesign():
    def filter_snps(self, df):
        """
        Iteratively check each SNP based on the specified rule and add it to a new DataFrame if it meets the condition:
        Keep a SNP if it is not less than 50 bp away from the previous SNP and not less than 300 bp away from the next SNP.

        :param df: pandas DataFrame containing the SNP data with columns 'Chr', 'position'
        :return: DataFrame with filtered SNPs.
        """
        # Sort the DataFrame by chromosome and position
        df = df.sort_values(['Chr', 'position'])

        # Calculate the distance to the next and previous SNP within each chromosome
        df['next_distance'] = df.groupby('Chr')['position'].shift(-1) - df['position']
        df['prev_distance'] = df['position'] - df.groupby('Chr')['position'].shift(1)

        # Set defaults for edge cases where there's no previous or next SNP
        df['prev_distance'] = df['prev_distance'].fillna(1000)  # No previous SNP, set previous distance to 1000
        df['next_distance'] = df['next_distance'].fillna(1000)  # No next SNP, set next distance to 1000

        # Create a new DataFrame to hold the filtered SNPs
        filtered_snps = []

        # Iterate through each row in the DataFrame
        for index, row in df.iterrows():
            # Apply the condition to decide whether to keep the SNP
            if (row['prev_distance'] > 50) and (row['next_distance'] > 300):
                # If the condition is met, add the row to the list
                filtered_snps.append(row)

        # Create a DataFrame from the list of kept rows
        filtered_df = pd.DataFrame(filtered_snps)

        return filtered_df

    def get_complement(self,sequence):
        complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join([complements[base] for base in sequence])  # Process each nucleotide

    def reverse_complement(self,seq):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement[base] for base in reversed(seq))



    def find_TM(self, seq):
        num_G = seq.count('G')
        num_C = seq.count('C')
        num_A = seq.count('A')
        num_T = seq.count('T')
        Tm = 64.9 + 41 * (num_G + num_C - 16.4) / (num_G + num_C + num_A + num_T)
        return Tm

    def find_reverse_primers(self, sequence, target_tm, tm_tolerance=2, min_length=20, max_length=30):
        potential_primers = []
        found = False  # Flag to indicate when two primers have been found
        for start in range(len(sequence)-30):
            for end in range(start + min_length, start + max_length+1):
                primer = sequence[start:end]
                reverse_tm = self.find_TM(primer)  # Corrected to pass 'primer' instead of 'seq'
                if target_tm - tm_tolerance <= reverse_tm <= target_tm + tm_tolerance:
                    potential_primers.append((primer, reverse_tm,200-end+1))
                    if len(potential_primers) >= 5:  # Stop if we already have 2 primers
                        found = True  # Set the flag to true as we've found enough primers
                        break  # Break out of the inner loop
            if found:  # Check the flag after breaking out of the inner loop
                break  # Break out of the outer loop if the flag is set
        return potential_primers[:1]  # Ensure to return only up to 1 primers, just in case

    def primer_dimer_check(self, primer1, primer2, threshold=5):
        # Check last 'threshold' bases for potential dimer formation
        end_primer1 = primer1[-threshold:]
        end_primer2 = self.reverse_complement(primer2[-threshold:])
        return end_primer1 == end_primer2
    def find_hairpins(self,seq, min_loop=3, min_stem=4):
        """Simple hairpin detection in a sequence."""
        for i in range(len(seq)):
            for j in range(i + min_stem * 2 + min_loop, len(seq) + 1):
                stem1 = seq[i:i + min_stem]
                #loop = seq[i + min_stem:i + min_stem + min_loop]
                stem2 = seq[j - min_stem:j]
                complement_stem2=self.get_complement(stem2)
                if stem1 == complement_stem2:
                    return True
        return False  # No hairpin structure found

    def calculate_gc_content(self, sequence):
        # Normalize the sequence to uppercase for consistent counting
        sequence = sequence.upper()
        # Count G and C bases
        gc_count = sequence.count('G') + sequence.count('C')
        # Calculate GC content as a percentage of the total length
        gc_content = (gc_count / len(sequence)) * 100
        return gc_content

    def find_repeats(self, sequence):
        """
        Find repeated motifs in the sequence based on predefined thresholds.
        Returns 'Yes' along with details of the exact number of repetitions if any significant repeats are found, otherwise 'No'.
        """
        # Predefined thresholds for each repeat length
        thresholds = {1: 6, 2: 4, 3: 4}  # Minimum repetitions: 6 for single, 5 for double, and 4 for triple nucleotides

        sequence_length = len(sequence)

        # Check for each repeat size up to max repeat length
        for n, min_repeats in thresholds.items():
            for i in range(sequence_length - n + 1):
                repeat = sequence[i:i + n]
                max_repeats = 0  # Variable to track the maximum exact number of repeats

                # Check if the next segment is the same as the current segment
                for j in range(1, (sequence_length - i) // n + 1):
                    if sequence[i:i + n * j] == repeat * j:
                        max_repeats = j  # Update the max_repeats if the current repeat count is valid
                    else:
                        break  # Break out of the loop once the sequence no longer matches the repeat pattern

                # After determining the maximum number of repeats, check if it meets or exceeds the threshold
                if max_repeats >= min_repeats:
                    return ("Yes",
                            f"{max_repeats} repetitions of a {n}-nucleotide repeat: '{repeat}'")  # Significant repeat found, return "Yes"

        return "No"  # No significant repeats found, return "No"
    def check(self,sequence,position,ref):
        SNP_ref=sequence[position]
        return SNP_ref==ref
    def KASP(self, sequences, chr, position, alter):
        # One-based position you're interested in
        target_position_one_based = position
        # Convert to zero-based for Python indexing
        target_position_zero_based = target_position_one_based - 1
        sequence = sequences.get(chr, None)
        check=self.check(sequence,target_position_zero_based,ref)
        for length in range(18, 35):  # Includes 18 to 30
            left_sequence = sequence[target_position_zero_based - length:target_position_one_based]
            Tm_forward = self.find_TM(left_sequence)
            # Check if Tm is within the desired range
            if 57 <= Tm_forward <= 65:
                # Construct primers with the selected left sequence
                A1 = FAM + left_sequence
                A2 = VIC + left_sequence[:-1] + alter
                #print(f"A1 primer: {A1}")
                #print(f"A2 primer: {A2}")

                # Proceed with the rest of your method
                right_sequence = sequence[target_position_one_based + 1:target_position_one_based + 200]
                #print(right_sequence)
                reversed_complementary_sequence = self.reverse_complement(right_sequence)
                potential_reverse_primers = self.find_reverse_primers(reversed_complementary_sequence,
                                                                      Tm_forward)
                results = []
                filtered_primers = []
                for primer, tm, pos in potential_reverse_primers:
                    filtered_primers.append((primer, tm, len(A1) + len(primer) + pos))
                    has_hairpin = self.find_hairpins(primer)
                    has_dimer = self.primer_dimer_check(A1, primer) or self.primer_dimer_check(A2, primer)
                    GC_forward = self.calculate_gc_content(A1)
                    GC_reverse = self.calculate_gc_content(primer)
                    find_repeat_F=self.find_repeats(A1)
                    find_repeat_R=self.find_repeats(primer)
                    if not has_hairpin and not has_dimer:
                        results.append((A1, A2, filtered_primers[-1:], "No", 'No',GC_forward,GC_reverse,find_repeat_F,find_repeat_R,check))
                    elif has_hairpin and not has_dimer:
                        results.append((A1, A2, filtered_primers[:1], "Yes", 'No',GC_forward,GC_reverse,find_repeat_F,find_repeat_R,check))
                    elif not has_hairpin and has_dimer:
                        results.append((A1, A2, filtered_primers[:1], "No", 'Yes',GC_forward,GC_reverse,find_repeat_F,find_repeat_R,check))
                    else:
                        results.append((A1, A2, filtered_primers[:1], "Yes", 'Yes',GC_forward,GC_reverse,find_repeat_F,find_repeat_R,check))
                    return results
                break



if __name__ == '__main__':
    # genome_path='C:\\genome\\2RBY\\Lcu.2RBY.FASTA'
    # md = MarkerDesign()  # Create an instance of the MarkerDesign class
    # marker_info = pd.read_csv('C:\genome\marker.csv')
    # output_csv_path = 'C:\\genome\\primer_designs.csv'
    # data = [] # Initialize a list to store data
    parser = argparse.ArgumentParser(description="Generate primers for KASP marker design.")
    parser.add_argument("-g", "--genome_path", required=True, help="Path to the genome FASTA file.")
    parser.add_argument("-m", "--marker_csv", required=True,
                        help="Path to the CSV file containing marker information.")
    parser.add_argument("-o", "--output_csv", required=True,
                        help="Path to save the output CSV file with primer designs.")
    args = parser.parse_args()
    # Read genome FASTA once and store sequences
    sequences = {}
    with open(args.genome_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            sequences[record.id] = record.seq

    md = MarkerDesign()
    marker_info = pd.read_csv(args.marker_csv)
    marker_filtered=md.filter_snps(marker_info)
    data = []

    for index, row in tqdm(marker_filtered.iterrows(), total=marker_filtered.shape[0], desc="Processing markers"):
    #for index, row in marker_filtered.iterrows():
        Chr = row.iloc[0]  # Selecting column 1
        position = row.iloc[1]  # Selecting column 2
        ref= row.iloc[2]
        alt = row.iloc[3]  # Selecting column 4
        result=md.KASP(sequences, Chr, position, alt)
        #print(result)
        if result is not None and result[0] is not None:
            A1_primer, A2_primer, reverse_primer, hairpin_check, dimer_check,GC_F,GC_R,repeat_F,repeat_R,check = result[0]
            primer, tm, size = reverse_primer[0]
            data.append({
                'Chr': Chr,
                'Position': position,
                'Ref':ref,
                'Alt': alt,
                'A1_Primer': A1_primer,
                'A2_Primer': A2_primer,
                'GC content for forward primer': GC_F,
                'Repeats within forward primer':repeat_F,
                'Reverse_Primer': primer,
                'GC content for reverse primer': GC_R,
                'Repeats within reverse primer': repeat_R,
                'Tm': tm,
                'Product size': size,
                'Hairpin':hairpin_check,
                'Primer dimer':dimer_check,
                'SNP alignment Check':check,
            })
        else:
            continue
    # Convert the collected data to a pandas DataFrame
    df = pd.DataFrame(data)

    # Save the DataFrame to a CSV file
    # Save the DataFrame to a CSV file
    df.to_csv(args.output_csv, index=False)
    print(f"Output saved to {args.output_csv}")
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
