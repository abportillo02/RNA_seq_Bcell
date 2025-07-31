#!/bin/bash
#SBATCH --job-name=d1_fastq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -p all
#SBATCH --mem=50G
#SBATCH --time=05:30:00






SRR_LIST=/home/abportillo/github_repo/RNA_seq_Bcell/scripts/SRR_Acc_List.txt

# Output directory (optional, create if needed)
OUTDIR=/home/abportillo/github_repo/RNA_seq_Bcell/scripts/fastq_bcell
mkdir -p "$OUTDIR"

module load SRA-Toolkit

# Loop through each SRR and run fasterq-dump
while read -r SRR; do
    echo "Processing $SRR"
    fasterq-dump "$SRR" -O "$OUTDIR"
done < "$SRR_LIST"