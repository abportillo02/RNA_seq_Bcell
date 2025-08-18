#!/bin/bash
#SBATCH --job-name=HERV_filter_gtf
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 16 # Number of cores
#SBATCH -N 1-4 # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --output=HERV_filter_gtf_%j.out

# source activate MyEnv  # if using conda/mamba
source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC


# Run the Python script
python /home/abportillo/github_repo/RNA_seq_Bcell/scripts/HERV_filter.py sampleNames.txt /home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/counts_tx /home/abportillo/github_repo/TEProf3/reference/journal.pcbi.1006453.s006