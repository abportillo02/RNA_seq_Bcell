#!/bin/sh

# In this round, we add up to 3 arguments in the config file, i.e. samples: sample_name path strand
# transcript level counts in this round3


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
  mkdir -p /home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/counts_tx
  outdir=/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/counts_tx
  ## make directory for each sample
  mkdir -p ${outdir}/${sample_name}
  ## software
  python=/home/abportillo/.conda/envs/mamba_abner_BC/bin/python
  samtools=/home/abportillo/.conda/envs/mamba_abner_BC/bin/samtools
  bamToBed=/home/abportillo/.conda/envs/mamba_abner_BC/bin/bamToBed
  coverageBed=/home/abportillo/.conda/envs/mamba_abner_BC/bin/coverageBed
  stringtie=/home/abportillo/.conda/envs/mamba_abner_BC/bin/stringtie
  prepDE_py=/home/abportillo/.conda/envs/mamba_abner_BC/bin/prepDE.py
  # annotations
  hg38_transcriptGTF=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/teAnno_round3/filtered_gencode_v46_chr_patch_hapl_scaff_annotation.gtf
  # write .sh script for each sample
  {
    echo -e "#!/bin/bash
#SBATCH --job-name=${sample_name}_getTxCounts
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 10 # Number of cores
#SBATCH -N 1-4 # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=30G
#SBATCH --time=06:00:00
#SBATCH --output=${outdir}/${sample_name}/${sample_name}_getTxCounts_%j.log\n

source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC


#### get number of aligned reads
${samtools} view -c -F 4 ${datapath_Bcell}/${sample_name}_sorted.bam \
> ${outdir}/${sample_name}/num_alignedReads.txt

#### transcript level counts ###########################################################################################

# #### get FPKM, TPM and read counts from stringtie outputs

${stringtie} ${datapath_Bcell}/${sample_name}_sorted.bam -e -B --rf \
-G ${hg38_transcriptGTF}  \
-o ${outdir}/${sample_name}/${sample_name}_Gencode_transcripts_ballgown.gtf -p 8

#### get TE initiated transcripts: FPKM, TPM

${python} /home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/hg38_p14/getTEtxReadCounts_fromGeneCodeTx.py \
${sample_name} \
${outdir} \
  /home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/hg38_p14/Gencode_TE_transcripts.bed

#### get Read Counts from stringtie outputs ############################################################################

echo '${sample_name} ${outdir}/${sample_name}/${sample_name}_Gencode_transcripts_ballgown.gtf' > \
${outdir}/${sample_name}/${sample_name}_ballgownGTFinfo.txt

##### use prepDE_py to get count matrix

${python} ${prepDE_py} -i ${outdir}/${sample_name}/${sample_name}_ballgownGTFinfo.txt \
-g ${outdir}/${sample_name}/${sample_name}_gene_count_matrix.csv \
-t ${outdir}/${sample_name}/${sample_name}_transcript_count_matrix.csv

##### get TE initiated transcripts: counts data

${python} /home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/hg38_p14/getTEtxReadCounts_fromGeneCodeTx.py \
${sample_name} \
${outdir} \
  /home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/hg38_p14/Gencode_TE_transcripts.bed
conda deactivate"

  } > ${outdir}/${sample_name}/${sample_name}_rnaGetTxCounts.sh
done < "${samples}"

echo "All sample script files created successfully."

