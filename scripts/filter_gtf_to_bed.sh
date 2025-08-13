#!/bin/bash
#SBATCH --job-name=filter_gtf
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 16 # Number of cores
#SBATCH -N 1-4 # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --output=filter_gtf_%j.out

# source activate MyEnv  # if using conda/mamba
source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC


# Run the Python script
python /home/abportillo/github_repo/RNA_seq_Bcell/scripts/getTEtx_fromGeneCodeTx.py
