import argparse
from Bio import SeqIO
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

FAM = 'GAAGGTGACCAAGTTCATGCT'
VIC = 'GAAGGTCGGAGTCAACGGATT'

class MarkerDesign():
    def get_complement(self, sequence):
        complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join([complements[base] for base in sequence])

    def reverse_complement(self, seq):
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
        found = False
        for start in range(len(sequence) - max_length):
            for end in range(start + min_length, start + max_length + 1):
                primer = sequence[start:end]
                reverse_tm = self.find_TM(primer)
                if target_tm - tm_tolerance <= reverse_tm <= target_tm + tm_tolerance:
                    potential_primers.append((primer, reverse_tm, 200 - end + 1))
                    if len(potential_primers) >= 2:
                        found = True
                        break
            if found:
                break
        return potential_primers[:1]
    def calculate_gc_content(self, sequence):
        # Normalize the sequence to uppercase for consistent counting
        sequence = sequence.upper()
        # Count G and C bases
        gc_count = sequence.count('G') + sequence.count('C')
        # Calculate GC content as a percentage of the total length
        gc_content = (gc_count / len(sequence)) * 100
        return gc_content
    def primer_dimer_check(self, primer1, primer2, threshold=5):
        end_primer1 = primer1[-threshold:]
        end_primer2 = self.reverse_complement(primer2[-threshold:])
        return end_primer1 == end_primer2

    def find_hairpins(self, seq, min_loop=3, min_stem=4):
        for i in range(len(seq)):
            for j in range(i + min_stem * 2 + min_loop, len(seq) + 1):
                stem1 = seq[i:i + min_stem]
                stem2 = seq[j - min_stem:j]
                complement_stem2 = self.get_complement(stem2)
                if stem1 == complement_stem2:
                    return True
        return False

    def KASP(self, sequences, chr, position, alter):
        target_position_one_based = position
        target_position_zero_based = target_position_one_based - 1
        sequence = sequences.get(chr, None)
        if sequence:
            for length in range(18, 35):
                left_sequence = sequence[target_position_zero_based - length:target_position_one_based]
                Tm_forward = self.find_TM(left_sequence)
                if 58 <= Tm_forward <= 63:
                    A1 = FAM + left_sequence
                    A2 = VIC + left_sequence[:-1] + alter
                    right_sequence = sequence[target_position_one_based + 1:target_position_one_based + 200]
                    reversed_complementary_sequence = self.reverse_complement(right_sequence)
                    potential_reverse_primers = self.find_reverse_primers(reversed_complementary_sequence, Tm_forward)
                    results = []
                    filtered_primers = []
                    for primer, tm, pos in potential_reverse_primers:
                        filtered_primers.append((primer, tm, len(A1) + len(primer) + pos))
                        has_hairpin = self.find_hairpins(primer)
                        has_dimer = self.primer_dimer_check(A1, primer) or self.primer_dimer_check(A2, primer)
                        GC_forward = self.calculate_gc_content(A1)
                        GC_reverse = self.calculate_gc_content(primer)
                        if not has_hairpin and not has_dimer:
                            results.append((A1, A2, filtered_primers[-1:], "No", 'No',GC_forward,GC_reverse))
                        elif has_hairpin and not has_dimer:
                            results.append((A1, A2, filtered_primers[:1], "Yes", 'No',GC_forward,GC_reverse))
                        elif not has_hairpin and has_dimer:
                            results.append((A1, A2, filtered_primers[:1], "No", 'Yes',GC_forward,GC_reverse))
                        else:
                            results.append((A1, A2, filtered_primers[:1], "Yes", 'Yes',GC_forward,GC_reverse))
                        return results
        return None

def process_marker(args_tuple):
    sequences, chr, position, alt, md = args_tuple
    return md.KASP(sequences, chr, position, alt)

def main():
    parser = argparse.ArgumentParser(description="Generate primers for KASP marker design.")
    parser.add_argument("-g", "--genome_path", required=True, help="Path to the genome FASTA file.")
    parser.add_argument("-m", "--marker_csv", required=True, help="Path to the CSV file containing marker information.")
    parser.add_argument("-o", "--output_csv", required=True, help="Path to save the output CSV file with primer designs.")
    args = parser.parse_args()

    # Read genome FASTA once and store sequences
    sequences = {}
    with open(args.genome_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            sequences[record.id] = record.seq

    md = MarkerDesign()
    marker_info = pd.read_csv(args.marker_csv)
    data = []
    task_args = [(sequences, row['Chr'], row['position'], row['alt'], md) for index, row in marker_info.iterrows()]

    with ThreadPoolExecutor() as executor:
        future_to_marker = {executor.submit(process_marker, arg): arg for arg in task_args}
        for future in tqdm(as_completed(future_to_marker), total=len(task_args), desc="Processing markers"):
            result = future.result()
            if result is not None and result[0] is not None:
                A1_primer, A2_primer, reverse_primer, hairpin_check, dimer_check,GC_F,GC_R = result[0]
                primer, tm, size = reverse_primer[0]
                data.append({
                    'Chr': future_to_marker[future][1],
                    'Position': future_to_marker[future][2],
                    'Alt': future_to_marker[future][3],
                    'A1_Primer': A1_primer,
                    'A2_Primer': A2_primer,
                    'GC content for forward primer': GC_F,
                    'Reverse_Primer': primer,
                    'GC content for reverse primer': GC_R,
                    'Tm': tm,
                    'Product size': size,
                    'Hairpin': hairpin_check,
                    'Primer dimer': dimer_check
                })

    df = pd.DataFrame(data)
    df.to_csv(args.output_csv, index=False)
    print(f"Output saved to {args.output_csv}")

if __name__ == '__main__':
    main()
