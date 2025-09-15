#!/bin/bash
#SBATCH --job-name=bed_to_bigwig_2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 16 # Number of cores
#SBATCH -N 1-4 # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --output=bed_to_bigwig_2.out

# Path to your directory with sorted BED files
dir="/home/abportillo/github_repo/RNA_seq_Bcell/output/chip_exo/KLF4/kzfp_peaks"

# Path to chromosome sizes file
chrom_sizes="/home/abportillo/github_repo/RNA_seq_Bcell/output/chip_exo/KLF4/kzfp_peaks/hg38.chrom.sizes"  

# Loop through each sorted BED file
for file in "${dir}"/*.sorted.bed; do
  base=$(basename "$file" .sorted.bed)
  bedGraphToBigWig "$file" "$chrom_sizes" "${dir}/${base}.bw"
done
