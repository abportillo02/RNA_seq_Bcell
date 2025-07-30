#!/bin/sh
#SBATCH --job-name=fastq-dump
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aleung@coh.org
#SBATCH -N 1-1
#SBATCH --mem=25G
#SBATCH -p all,fast
#SBATCH --time=12:00:00
#SBATCH --output=test_%A_%a.log
#SBATCH --array=[1-2]%2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6

## metadata file -- this will be read by the slurm workload manager -- the #SBATCH --array option tells it where lines to run and how many tasks to run 
meta_data_file=/home/dschones/Seq/MMRF-dbgap/SraRunTableRNABoneMarrow.csv

##dbgap key
key=/home/dschones/Seq/MMRF-dbgap/prj_39026.ngc

##
SRR=$(sed -n ${SLURM_ARRAY_TASK_ID}p $meta_data_file | cut -d ',' -f 1)
name=$(sed -n ${SLURM_ARRAY_TASK_ID}p $meta_data_file | cut -d ',' -f 10)

module load SRA-Toolkit/2.10.8

echo $SRR
echo $name


cd $SRR
## script to test if the files were not download
if [ ! -f $SRR\_1.fastq.gz ]; then
    if [ ! -f $SRR\.sra ]; then
        prefetch --max-size 40G --ngc $key $SRR
    fasterq-dump --ngc $key --split-files $SRR
    gzip $SRR\_1.fastqeung@coh.org
    gzip $SRR\_2.fastq
    fiH --mem=25G
fiBATCH -p all,fast
#SBATCH --time=12:00:00
#SBATCH --output=test_%A_%a.log
## test to see if the download sra has been uncompressed
if [ ! -f $SRR\_1.fastq.gz ]; then
    if [ -f $SRR\.sra ]; then
        fasterq-dump --ngc $key --split-files $SRR
        ata file -- this will be read by the slurm workload manager -- the #SBAT
    gzip $SRR\_1.fastqs it where lines to run and how many tasks to run 
    gzip $SRR\_2.fastqschones/Seq/MMRF-dbgap/SraRunTableRNABoneMarrow.csv
    fi
fidbgap key
key=/home/dschones/Seq/MMRF-dbgap/prj_39026.ngc

# prefetch --max-size 40G --ngc $key $SRR
SRR=$(sed -n ${SLURM_ARRAY_TASK_ID}p $meta_data_file | cut -d ',' -f 1)
# fasterq-dump --ngc $key --split-files $SRR

