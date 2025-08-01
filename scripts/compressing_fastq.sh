#!/bin/sh
#SBATCH --job-name=fastq_compress
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -p all
#SBATCH --mem=50G
#SBATCH --time=05:30:00
#SBATCH --output=fastq_compress.log


# Directory containing the files
fastq_compress="/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell"
# Compress all .fastq files in parallel and keep the originals
find "$fastq_compress" -type f -name "*.fastq" -print0 | xargs -0 -P 4 -n 1 gzip -N
