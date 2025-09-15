#!/bin/bash
#SBATCH --job-name=homer_motif_discovery
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -p all
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --output=/home/abportillo/github_repo/RNA_seq_Bcell/output/chip_exo/KLF4/homer_output/homer_motif.log

# === User-defined variables ===
FASTA_FILE="/home/abportillo/github_repo/RNA_seq_Bcell/output/chip_exo/KLF4/KZFP519_overlapping_KLF4.fa"
OUTPUT_DIR="/home/abportillo/github_repo/RNA_seq_Bcell/output/chip_exo/KLF4/homer_output"
MOTIF_LENGTHS="8,10,12,14"
NUM_THREADS=4

source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC
module load homer

echo "Running HOMER motif discovery on $FASTA_FILE..."
findMotifs.pl "$FASTA_FILE" fasta "$OUTPUT_DIR" -len "$MOTIF_LENGTHS" -p "$NUM_THREADS"
echo "Motif discovery complete. Results saved in $OUTPUT_DIR"
