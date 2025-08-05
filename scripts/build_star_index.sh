#!/bin/bash

#SBATCH --job-name=Build_STAR_Index
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 16 # Number of cores
#SBATCH -N 1-4 # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --output=Build_STAR_Index.log


cd /home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/
mkdir -p STAR_hg38_p14_geneCodeGTF_filter

# get filtered geneCode annotation, filter out patches
zcat gencode_v46_chr_patch_hapl_scaff_annotation.gtf.gz | \
awk '$1 ~ /^chr([1-9]$|1[0-9]$|2[0-2]$|X$|Y$|M$)/' > filtered_gencode_v46_chr_patch_hapl_scaff_annotation.gtf

# check the filtered geneCode annotation
cut -f1 filtered_gencode_v46_chr_patch_hapl_scaff_annotation.gtf | sort | uniq

cd STAR_hg38_p14_geneCodeGTF_filter
# build STAR index
/home/abportillo/.conda/envs/mamba_abner_BC/bin/STAR \
--runMode genomeGenerate \
--runThreadN 8 \
--genomeDir /home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/hg38_p14/STAR_hg38_p14_geneCodeGTF_filter \
--genomeFastaFiles /home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/hg38_p14/filtered_hg38_p14.fa \
--sjdbGTFfile /home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/hg38_p14/filtered_gencode_v46_chr_patch_hapl_scaff_annotation.gtf \
--sjdbOverhang 100

ln -s /home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/hg38_p14/filtered_gencode_v46_chr_patch_hapl_scaff_annotation.gtf .