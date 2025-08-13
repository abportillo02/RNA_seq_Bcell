#!/bin/bash
#SBATCH --job-name=teprof3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 36 # Number of cores
#SBATCH -N 1-4 # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=teprof3.log

# R1: use --quanreadlength 75

source /home/abportillo/.bashrc
module load Mamba/24.3.0-0
mamba activate /home/abportillo/.conda/envs/mamba_abner_BC

# Define the directory:
teprof3Dir="/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/TEProf3_TEdetect"
cd ${teprof3Dir}

# Clean old links and manifest
# rm -f *.bam *.bam.bai *.tab *.gtf sample_manifest.txt



# Create symbolic links to the necessary files
ln -s /home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/*/*_sorted_nr_sorted.bam .
ln -s /home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/*/*_sorted_nr_sorted.bam.bai .
ln -s /home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/*/*_SJ.out.tab .
ln -s /home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/counts_tx/*/*_Gencode_transcripts_ballgown.gtf .

# get sample_manifest file
rm sample_manifest.txt
find . -maxdepth 1 -name "*.bam"  | while read file ; do xbase=$(basename $file); sample_name=${xbase/_sorted_nr_sorted.bam/} ; echo -e "${sample_name}\tshort\t${xbase}" >> sample_manifest.txt; done ;
find . -maxdepth 1 -name "*.tab"  | while read file ; do xbase=$(basename $file); sample_name=${xbase/_SJ.out.tab/} ; echo -e "${sample_name}\tSJ\t${xbase}" >> sample_manifest.txt; done ;
find . -maxdepth 1 -name "*.gtf"  | while read file ; do xbase=$(basename $file); sample_name=${xbase/.gtf/} ; echo -e "${sample_name}\tgtf\t${xbase}" >> sample_manifest.txt; done

# reset all folders if error log
teprof3 -f sample_manifest.txt --reset

# run teprof3:
# Assemble step (-as) \ # run Process assemble step (-ps) \ # run Filter transcripts (-fs) \\ 
# run Mega assembly (--tacothread) \ # run Quantification (-qs)
teprof3 -f sample_manifest.txt -s 11 -am 1 -as 11 --assemblethread 4 --assemblestrand 1 --assemblejunctionread 2 \
  -ps 11 -pt 0.5 -ptn 100 \
  -fs 11 -fm 1 --filterintronretention 3 \
  --filtermonoexontpm 1 --filterdownstreammate 2 --filterratio 0.5 \
  --tacothread 11 \
  -qs 11 --quansamplenumbercon 100 --quanmode 1 --quanreadlength 75

mamba deactivate

