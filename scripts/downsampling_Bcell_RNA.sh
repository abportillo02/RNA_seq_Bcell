#!/bin/bash

# use `seqtk` software to downsampling PMM fastq files to 20000 reads
module load seqtk

# datapath
datapath_Bcell=/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell

# output dir
outdir=/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/downsample_Bcell_RNA
cd ${outdir}

paste samples.txt <(seq 100 110) > seeded_samples.txt

# get 20000 reads from each sample (PMM)
while read -r sample_name seed; do
  mkdir -p ${sample_name}

  seqtk sample -s ${seed} ${datapath_Bcell}/${sample_name}*_R1*.fastq.gz 20000 | gzip > \
  ${outdir}/${sample_name}/${sample_name}_R1_sub20000.fq.gz

  seqtk sample -s ${seed} ${datapath_Bcell}/${sample_name}*_R2*.fastq.gz 20000 | gzip > \
  ${outdir}/${sample_name}/${sample_name}_R2_sub20000.fq.gz

done < seeded_samples.txt 
