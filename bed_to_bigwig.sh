#!/bin/bash
#SBATCH --job-name=bed_to_bigwig
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 16 # Number of cores
#SBATCH -N 1-4 # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --output=bed_to_bigwig.out

for file in /home/abportillo/github_repo/RNA_seq_Bcell/scripts/chip_exo/KLF4/kzfp_peaks/*.bed; do
  base=$(basename "$file" .bed)
  sort -k1,1 -k2,2n "$file" > "/home/abportillo/github_repo/RNA_seq_Bcell/scripts/chip_exo/KLF4/kzfp_peaks/${base}.sorted.bed"
  bedGraphToBigWig "/home/abportillo/github_repo/RNA_seq_Bcell/scripts/chip_exo/KLF4/kzfp_peaks/${base}.sorted.bed" hg38.chrom.sizes "/home/abportillo/github_repo/RNA_seq_Bcell/scripts/chip_exo/KLF4/kzfp_peaks/${base}.bw"
done
