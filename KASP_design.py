import argparse
from Bio import SeqIO
import pandas as pd
from tqdm import tqdm  # Import tqdm for the progress bar
FAM = 'GAAGGTGACCAAGTTCATGCT'
VIC = 'GAAGGTCGGAGTCAACGGATT'
class MarkerDesign():

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
        for start in range(len(sequence)-30):
            if len(potential_primers) >= 5:  # Stop if we already have 5 primers
                break
            for end in range(start + min_length, start + max_length+1):

                primer = sequence[start:end]
                reverse_tm = self.find_TM(primer)  # Corrected to pass 'primer' instead of 'seq'
                if target_tm - tm_tolerance <= reverse_tm <= target_tm + tm_tolerance:
                    potential_primers.append((primer, reverse_tm,start))
        return potential_primers[:5]  # Ensure to return only up to 5 primers, just in case

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
                loop = seq[i + min_stem:i + min_stem + min_loop]
                stem2 = seq[j - min_stem:j]
                complement_stem2=self.get_complement(stem2)
                if stem1 == complement_stem2:
                    return True
        return False  # No hairpin structure found
    def KASP(self, genome_path, chr, position, alter):
        # One-based position you're interested in
        target_position_one_based = position
        # Convert to zero-based for Python indexing
        target_position_zero_based = target_position_one_based - 1
        for record in SeqIO.parse(genome_path, 'fasta'):
            if record.id == chr:
                # Iterate through possible left sequence lengths (18 to 30)
                for length in range(18, 35):  # Includes 18 to 30
                    left_sequence = record.seq[target_position_zero_based - length:target_position_one_based]
                    Tm_forward = self.find_TM(left_sequence)
                    # Check if Tm is within the desired range
                    if 58 <= Tm_forward <= 63:
                        # Construct primers with the selected left sequence
                        A1 = FAM + left_sequence
                        A2 = VIC + left_sequence[:-1] + alter
                        #print(f"A1 primer: {A1}")
                        #print(f"A2 primer: {A2}")

                        # Proceed with the rest of your method
                        right_sequence = record.seq[target_position_one_based + 1:target_position_one_based + 200]
                        reversed_complementary_sequence = self.reverse_complement(right_sequence)
                        potential_reverse_primers = self.find_reverse_primers(reversed_complementary_sequence,
                                                                              Tm_forward)
                        results = []
                        filtered_primers = []
                        for primer, tm, pos in potential_reverse_primers:
                            filtered_primers.append((primer, tm, len(A1) + len(A2) + pos))
                            has_hairpin = self.find_hairpins(primer)
                            has_dimer = self.primer_dimer_check(A1, primer) or self.primer_dimer_check(A2, primer)
                            if not has_hairpin and not has_dimer:
                                results.append((A1, A2, filtered_primers[-1:], "No", 'No'))
                            elif has_hairpin and not has_dimer:
                                results.append((A1, A2, filtered_primers[:1], "Yes", 'No'))
                            elif not has_hairpin and has_dimer:
                                results.append((A1, A2, filtered_primers[:1], "No", 'Yes'))
                            else:
                                results.append((A1, A2, filtered_primers[:1], "Yes", 'Yes'))
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

    md = MarkerDesign()
    marker_info = pd.read_csv(args.marker_csv)
    data = []
    for index, row in tqdm(marker_info.iterrows(), total=marker_info.shape[0], desc="Processing markers"):
    #for index, row in marker_info.iterrows():
        Chr = row.iloc[0]  # Selecting column 1
        position = row.iloc[1]  # Selecting column 2
        alt = row.iloc[2]  # Selecting column 3
        A1_primer, A2_primer, reverse_primer, hairpin_check, dimer_check = md.KASP(args.genome_path, Chr, position, alt)[0]
        primer, tm, size = reverse_primer[0]
        data.append({
            'Chr': Chr,
            'Position': position,
            'Alt': alt,
            'A1_Primer': A1_primer,
            'A2_Primer': A2_primer,
            'Reverse_Primer': primer,
            'Tm': tm,
            'Product size': size,
            'hairpin':hairpin_check,
            'Primer dimer':dimer_check
        })
    # Convert the collected data to a pandas DataFrame
    df = pd.DataFrame(data)

    # Save the DataFrame to a CSV file
    # Save the DataFrame to a CSV file
    df.to_csv(args.output_csv, index=False)
    print(f"Output saved to {args.output_csv}")
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
