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

# read each line from the $1 input sampleNames.txt file and create .sh script for each sample
while IFS=" " read -r sample_name path; do
  echo "Creating bash script for sample: $sample_name"
  ## data path
  datapath_Bcell=/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/${sample_name}
  mkdir -p /home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/gene_counts
  outdir=/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/gene_counts
  ## make directory for each sample
  mkdir -p ${outdir}/${sample_name}
  ## software
  featureCounts=/home/abportillo/.conda/envs/mamba_abner_BC/bin/featureCounts
  # annotations
  hg38_transcriptGTF=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/teAnno_round3/filtered_gencode_v46_chr_patch_hapl_scaff_annotation.gtf
#  hg38_mRNAexonsSAF=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14_saf/hg38_p14_mRNA.saf
  # write .sh script for each sample
  {
    echo -e "#!/bin/bash
#SBATCH --job-name=${sample_name}_getGeneCounts
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 10 # Number of cores
#SBATCH -N 1-4 # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=30G
#SBATCH --time=06:00:00
#SBATCH --output=${outdir}/${sample_name}_getGeneCounts_%j.log\n

source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC

#### gene level counts ###########################################################################################

${featureCounts} -B -p --countReadPairs -Q 30 -s 2 -T 8 -a ${hg38_transcriptGTF} \
-o ${outdir}/${sample_name}/${sample_name}_gene_counts.txt \
${datapath_Bcell}/${sample_name}_crick_merged.bam ${datapath_Bcell}/${sample_name}_watson_merged.bam

conda deactivate"
  } > ${outdir}/${sample_name}/${sample_name}_rnaGetGeneCounts.sh
#  cd ${outdir}/${sample_name}
#  sbatch ${sample_name}_rnaGetGeneCounts.sh
done < "${samples}"

echo "All sample script files created successfully."
