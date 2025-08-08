#!/bin/bash


if [ -z "$1" ]; then
  echo "Usage: $0 <sampleNames.txt>"
  exit 1
fi

sample_list=$1
outdir="/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/Transcript_finder"

mkdir -p "${outdir}"

while IFS= read -r sample; do
  [ -z "$sample" ] && continue

  echo "Creating bash script for sample: $sample"

  sample_dir="/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/${sample}"

  {
  echo -e "#!/bin/bash
#SBATCH --job-name=${sample}_TEProf3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 16
#SBATCH -N 1-4
#SBATCH -p all
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --output=${outdir}/${sample}_TEProf3.log

cd ${sample_dir}

source /home/abportillo/.bashrc
module load Mamba/24.3.0-0
mamba activate /home/abportillo/.conda/envs/mamba_abner_BC

# rm -f sample_manifest.txt
find . -maxdepth 1 -name "*.bam" | while read file; do xbase=\$(basename \$file); sample_name=\${xbase/_sorted_nr_sorted.bam/}; echo -e \"\${sample_name}\tshort\t\${xbase}\" >> sample_manifest.txt; done
find . -maxdepth 1 -name "*.tab" | while read file; do xbase=\$(basename \$file); sample_name=\${xbase/_SJ.out.tab/}; echo -e \"\${sample_name}\tSJ\t\${xbase}\" >> sample_manifest.txt; done
find . -maxdepth 1 -name "*.gtf" | while read file; do xbase=\$(basename \$file); sample_name=\${xbase/.gtf/}; echo -e \"\${sample_name}\tgtf\t\${xbase}\" >> sample_manifest.txt; done

# Uncomment if you want to reset before running
# teprof3 -f sample_manifest.txt --reset

teprof3 -f sample_manifest.txt -s 10 -am 1 -as 10 --assemblethread 4 --assemblestrand 1 --assemblejunctionread 2 \\
  -ps 10 -pt 0.5 -ptn 100 \\
  -fs 10 -fm 1 --filterintronretention 3 \\
  --filtermonoexontpm 1 --filterdownstreammate 2 --filterratio 0.5 \\
  --tacothread 10 \\
  -qs 10 --quansamplenumbercon 100 --quanmode 1 --quanreadlength 75

mamba deactivate
"
  } > "${outdir}/${sample}_runTEProf3.sh"

  echo "Created script: ${outdir}/${sample}_runTEProf3.sh"

done < "$sample_list"

echo "Finished creating scripts for all samples in $sample_list"
