#!/bin/bash
#SBATCH --job-name=bed_to_bigwig
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 16
#SBATCH -N 1-4
#SBATCH -p all
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --output=bed_to_bigwig.out

# Directory containing BED files
dir="/home/abportillo/github_repo/RNA_seq_Bcell/output/chip_exo/KLF4/kzfp_peaks"

# Path to chromosome sizes file
chrom_sizes="/home/abportillo/genomes/hg38.chrom.sizes"  # Update if needed

# Loop through each BED file
for file in "${dir}"/*.bed; do
  base=$(basename "$file" .bed)

  # Extract only the required columns: chr, start, end, score
  awk '{print $1"\t"$2"\t"$3"\t"$5}' "$file" > "${dir}/${base}.bedgraph"

  # Sort the bedGraph file
  sort -k1,1 -k2,2n "${dir}/${base}.bedgraph" > "${dir}/${base}.sorted.bedgraph"

  # Convert to BigWig
  bedGraphToBigWig "${dir}/${base}.sorted.bedgraph" "$chrom_sizes" "${dir}/${base}.bw"
done
