import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class ChromosomeSNP():
    def load_snp(self, snp_file):
        df = pd.read_csv(snp_file)
        df = df.sort_values(['Chr', 'Position'])
        return df

    def draw_chromosomes_and_snps(self, chromosome_lengths, file_path):
        snp_info = self.load_snp(file_path)
        chromosome_lengths = np.array(chromosome_lengths) / 1000000
        # Determine the number of chromosomes
        n = len(chromosome_lengths)

        # Create a figure and an axis for plotting
        fig, ax = plt.subplots(figsize=(10, 8))

        # Positioning each chromosome on the plot
        y_positions = np.linspace(3, 10 * n, n)

        # Drawing chromosomes and SNPs
        for i, length in enumerate(chromosome_lengths):
            # Draw the chromosome as a horizontal line
            x_start = 0
            x_end = x_start + length
            print(x_end)
            y = y_positions[i]
            ax.plot([x_start, x_end], [y, y], color='gray', alpha=0.7, linewidth=8)

            # Filter SNP data for this chromosome
            chr_label = f"Lcu.2RBY.Chr{i + 1}"
            chr_snps = snp_info[snp_info['Chr'] == chr_label]
            snp_positions = chr_snps['Position'].values/1000000
            # Draw vertical lines for SNPs
            for snp_position in snp_positions:
                ax.vlines(x=snp_position, ymin=y - 1, ymax=y + 1, colors='red', linewidth=0.3)
            # Draw SNPs as red points
            #ax.scatter(snp_positions, [y] * len(snp_positions), color='red', s=7, label='SNPs' if i == 0 else "")

        # Set the limits and labels
        ax.set_xlim(0, max(chromosome_lengths) + 20)
        ax.set_ylim(0, 10 * (n + 1))
        ax.set_xlabel("Position (MB)")
        ax.set_ylabel("Chromosome")
        ax.set_yticks(y_positions)
        ax.set_yticklabels([f"Lcu.2RBY.Chr{i + 1}" for i in range(n)])
        #ax.legend()

        # Hide the axes for a cleaner look
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        plt.savefig('C:\genome\LR-66_marker_design\SNP_density.jpg', format='jpeg', dpi=600)
        plt.show()

# Example usage
chromosome_lengths = [538362633, 614089212, 430115100, 482216434, 475007465, 420584759, 529090306]  # Example lengths for chromosomes
file_path = 'C:\genome\LR-66_marker_design\primers_LR-66-590.csv'  # Replace this with your actual file path

# Create an instance of the class and call the plotting function
chromosome_plotter = ChromosomeSNP()
chromosome_plotter.draw_chromosomes_and_snps(chromosome_lengths, file_path)
