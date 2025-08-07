#!/bin/sh

# run allo after deduplication in this round of analysis, and consider SJ (splice introns out when
# summing uniquely mapped reads)
# use TCGA 2passSTAR alignment and use hg38.p14 without patches

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
  
  # software
  # fastqc=/home/qwan/miniconda3/envs/coh/bin/fastqc
  # fastp=/home/qwan/miniconda3/envs/coh/bin/fastp
  # allo=/home/qwan/miniconda3/envs/coh/bin/allo
  # samtools=/home/qwan/miniconda3/envs/coh/bin/samtools
  # java=/home/qwan/miniconda3/envs/coh/bin/java
  # bedtools=/home/qwan/miniconda3/envs/coh/bin/bedtools

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
### Align reads with STAR:
${STAR} --genomeDir /home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/hg38_p14/STAR_hg38_p14_geneCodeGTF_filter \
--readFilesIn ${datapath_Bcell}/${sample_name}_1.fastq.gz ${datapath_Bcell}/${sample_name}_2.fastq.gz \
--readFilesCommand zcat \
--runThreadN 8 \
--twopassMode Basic \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outFilterType BySJout \
--outFilterScoreMinOverLread 0.33 \
--outFilterMatchNminOverLread 0.33 \
--limitSjdbInsertNsj 1200000 \
--outFileNamePrefix ${outdir}/${sample_name}_ \
--outSAMstrandField intronMotif \
--outFilterIntronMotifs None \
--alignSoftClipAtReferenceEnds Yes \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within \
--genomeLoad NoSharedMemory \
--chimSegmentMin 15 \
--chimJunctionOverhangMin 15 \
--chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
--chimOutJunctionFormat 1 \
--chimMainSegmentMultNmax 1 \
--outSAMattributes NH HI AS nM NM ch

#--outSAMattrRGline ID:sample1 SM:sample1 PL:ILLUMINA LB:lib1 PU:unit1; not sure lib info. so leave it out
#--outFilterMismatchNmax 10, not use in R1 and R2 because:
# other para. covered this filtering, the default value is unlimited
#--outFilterMultimapNmax 100, for repeat, not use in this round
#--winAnchorMultimapNmax 100, for repeat, not use in this round
# --outSAMmultNmax 25, for allo, not use in this round

### sort by coordinates, Remove duplicates and sort

${samtools} sort -@ 8 -O bam -o ${outdir}/${sample_name}_sorted.bam ${outdir}/${sample_name}_Aligned.out.bam
${samtools} index ${outdir}/${sample_name}_sorted.bam
rm ${outdir}/${sample_name}_Aligned.out.bam

### sort by coordinates, Remove duplicates and sort
# module load picard/2.21.1
picard MarkDuplicates \\
    I=${outdir}/${sample_name}_sorted.bam \\
    O=${outdir}/${sample_name}_nr_sorted.bam  \\
    M=${outdir}/${sample_name}_picardStats.txt \\
    REMOVE_DUPLICATES=true 
   
picard AddOrReplaceReadGroups \
    -I ${outdir}/${sample_name}_sorted.bam \
    -O ${outdir}/${sample_name}_rg_sorted.bam \
    --RGID ${sample_name} \
    --RGLB default_lib \
    --RGPL ILLUMINA \
    --RGPU unknown \
    --RGSM ${sample_name}\
    --SORT_ORDER coordinate

 ${java} -Djava.io.tmpdir=/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/temp \
  -jar /home/abportillo/.conda/envs/mamba_abner_BC/share/picard-3.3.0-0/picard.jar MarkDuplicates \
  --INPUT ${outdir}/${sample_name}_rg_sorted.bam --OUTPUT ${outdir}/${sample_name}_nr_sorted.bam \
 --REMOVE_DUPLICATES true --READ_NAME_REGEX null --METRICS_FILE ${outdir}/${sample_name}_picardStats.txt

${samtools} sort -@ 8 -O bam -o ${outdir}/${sample_name}_sorted_nr_sorted.bam ${outdir}/${sample_name}_nr_sorted.bam
${samtools} index -@ 8 ${outdir}/${sample_name}_sorted_nr_sorted.bam

### check for strand info. of RNAseq data here, if stranded, get separate bam for watson and crick
module load RSeQC/4.0.0
cd ${outdir}
infer_experiment.py \
-r /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/teAnno_round3/hg38_gencode_exon.bed \
-i ${outdir}/${sample_name}_sorted_nr_sorted.bam

### get bam files on watson and crick strands
${samtools} view -@ 7 -h -f 99 -S -b ${outdir}/${sample_name}_sorted_nr_sorted.bam > ${outdir}/crick1.bam
${samtools} view -@ 7 -h -f 147 -S -b ${outdir}/${sample_name}_sorted_nr_sorted.bam > ${outdir}/crick2.bam
${samtools} view -@ 7 -h -f 83 -S -b ${outdir}/${sample_name}_sorted_nr_sorted.bam > ${outdir}/watson1.bam
${samtools} view -@ 7 -h -f 163 -S -b ${outdir}/${sample_name}_sorted_nr_sorted.bam > ${outdir}/watson2.bam
${samtools} sort -@ 7 -o ${outdir}/${sample_name}_crick1_sorted.bam ${outdir}/crick1.bam
${samtools} sort -@ 7 -o ${outdir}/${sample_name}_crick2_sorted.bam ${outdir}/crick2.bam
${samtools} sort -@ 7 -o ${outdir}/${sample_name}_watson1_sorted.bam ${outdir}/watson1.bam
${samtools} sort -@ 7 -o ${outdir}/${sample_name}_watson2_sorted.bam ${outdir}/watson2.bam

${samtools} merge -f ${outdir}/${sample_name}_crick_merged.bam \
${outdir}/${sample_name}_crick1_sorted.bam ${outdir}/${sample_name}_crick2_sorted.bam
${samtools} merge -f ${outdir}/${sample_name}_watson_merged.bam \
${outdir}/${sample_name}_watson1_sorted.bam ${outdir}/${sample_name}_watson2_sorted.bam

${samtools} index ${outdir}/${sample_name}_crick_merged.bam
${samtools} index ${outdir}/${sample_name}_watson_merged.bam

rm ${outdir}/crick1.bam ${outdir}/crick2.bam ${outdir}/watson1.bam ${outdir}/watson2.bam
rm ${outdir}/${sample_name}_crick1_sorted.bam ${outdir}/${sample_name}_crick2_sorted.bam \
${outdir}/${sample_name}_watson1_sorted.bam ${outdir}/${sample_name}_watson2_sorted.bam

### bam to bedGraph
${bedtools} genomecov -bg -ibam ${outdir}/${sample_name}_crick_merged.bam \
-g /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/filtered_hg38_p14.chrom.sizes \
> ${outdir}/${sample_name}_crick_readCount.bedGraph

${bedtools} genomecov -bg -ibam ${outdir}/${sample_name}_watson_merged.bam \
-g /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/filtered_hg38_p14.chrom.sizes \
> ${outdir}/${sample_name}_watson_readCount.bedGraph

### bedGraph to bigwig (read counts)
${wigToBigWig} ${outdir}/${sample_name}_crick_readCount.bedGraph \
/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/filtered_hg38_p14.chrom.sizes \
${outdir}/${sample_name}_crick_readCount.bw

${wigToBigWig} ${outdir}/${sample_name}_watson_readCount.bedGraph \
/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/filtered_hg38_p14.chrom.sizes \
${outdir}/${sample_name}_watson_readCount.bw

### bigwig (cpm)
${bamCoverage} -b ${outdir}/${sample_name}_crick_merged.bam \
-o ${outdir}/${sample_name}_crick_cpm.bw --normalizeUsing CPM \
--numberOfProcessors 12 \
--effectiveGenomeSize 2913022398 --binSize 1

${bamCoverage} -b ${outdir}/${sample_name}_watson_merged.bam \
-o ${outdir}/${sample_name}_watson_cpm.bw --normalizeUsing CPM \
--numberOfProcessors 12 \
--effectiveGenomeSize 2913022398 --binSize 1

### bigwig (rpkm)
${bamCoverage} -b ${outdir}/${sample_name}_crick_merged.bam \
-o ${outdir}/${sample_name}_crick_rpkm.bw --normalizeUsing RPKM \
--numberOfProcessors 12 \
--effectiveGenomeSize 2913022398 --binSize 1

${bamCoverage} -b ${outdir}/${sample_name}_watson_merged.bam \
-o ${outdir}/${sample_name}_watson_rpkm.bw --normalizeUsing RPKM \
--numberOfProcessors 12 \
--effectiveGenomeSize 2913022398 --binSize 1

rm ${outdir}/*.bedGraph
conda deactivate"
  } > ${outdir}/${sample_name}_rnaPreprocess.sh
#  cd ${outdir}
#  sbatch ${sample_name}_rnaPreprocess.sh
done < ${samples}

echo "All sample_name script files created successfully."





