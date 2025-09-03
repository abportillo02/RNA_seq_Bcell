#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --output=sra_download.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 16 # Number of cores
#SBATCH -N 1-4 # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=150G
#SBATCH --time=48:00:00

# Load Conda environment
source source /home/abportillo/.bashrc
conda activate mamba_abner_BC

# Output directory
output_dir="/home/abportillo/github_repo/RNA_seq_Bcell/scripts/methylation"
mkdir -p "$output_dir"

# SRA IDs
sra_ids=(
  SRR948841 SRR948842 SRR948843 SRR948844 SRR948845 SRR948846
  SRR948847 SRR948848 SRR948849 SRR948850 SRR948851 SRR948852
  SRR948853 SRR948854 SRR948855
)

# Loop through and download each
for sra in "${sra_ids[@]}"; do
  echo "Downloading $sra..."
  fasterq-dump "$sra" -O "$output_dir" -e 4
done

echo "All downloads complete."
