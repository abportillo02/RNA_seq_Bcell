
#!/bin/sh

#Fastqc the fastq files 
# Then map data to hg38 and get ERV loci 

# Check if the input file is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <input_file>"
  exit 1
fi

## variables
samples=$1

# Check if the input file exists
if [ ! -f "$samples" ]; then
  echo "Input file not found!"
  exit 1
fi


# read each line from the $1 input sampleNames.txt file and create .sh script for each sample_name
while IFS= read -r sample_name; do
  echo "Creating bash script for sample_name: $sample_name"
  ## data path
  ## we can soft-link all data to a data folder, so we can use it more conveniently
  datapath_Bcell=/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell
  mkdir -p /home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess
  mkdir -p /home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/${sample_name}
  outdir=/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/${sample_name}

  #software paths
  java=/home/abportillo/.conda/envs/mamba_abner_BC/bin/java
  bamCoverage=/home/abportillo/.conda/envs/mamba_abner_BC/bin/bamCoverage
  samtools=/home/abportillo/.conda/envs/mamba_abner_BC/bin/samtools
  STAR=/home/abportillo/.conda/envs/mamba_abner_BC/bin/STAR
  picard=/home/abportillo/.conda/envs/mamba_abner_BC/bin/picard
  wigToBigWig=/home/abportillo/.conda/envs/mamba_abner_BC/bin/wigToBigWig

  {
    echo -e "#!/bin/bash
#SBATCH --job-name=${sample_name}_RNA_hg38_p14_2passStar
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 16 # Number of cores
#SBATCH -N 1-4 # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --output=${outdir}/${sample_name}_RNA_hg38_p14_2passStar_%j.log

source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC

#### fastqc for fastq files
mkdir -p ${outdir}/fastqc_out
module load FastQC/0.11.8
fastqc -t 8 -o ${outdir}/fastqc_out \
${datapath_Bcell}/${sample_name}_1.fastq.gz ${datapath_Bcell}/${sample_name}_2.fastq.gz

ln -s ${datapath_Bcell}/${sample_name}_1.fastq.gz ${outdir}/${sample_name}_R1.fastq.gz
ln -s ${datapath_Bcell}/${sample_name}_2.fastq.gz ${outdir}/${sample_name}_R2.fastq.gz
module unload FastQC/0.11.8


#### Mapping to hg38/ERVmap


## Setup for pair-end 
