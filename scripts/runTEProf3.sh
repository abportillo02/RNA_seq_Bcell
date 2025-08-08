#!/bin/bash

# Check if file with folder names is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <sampleNames.txt>"
  exit 1
fi

folder_list=$1
outdir="/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/Transcript_finder"  # single shared folder

# Make sure the one output directory exists
mkdir -p "${outdir}"

while IFS= read -r folder; do
  # Skip empty lines
  [ -z "$folder" ] && continue

  echo "Creating bash script for folder: $folder"
  
  {
  echo -e "#!/bin/bash
#SBATCH --job-name=${folder}_TEProf3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 16 # Number of cores
#SBATCH -N 1-4 # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --output=${outdir}/${folder}_TEProf3.log

cd ${outdir} 

source /home/abportillo/.bashrc
module load Mamba/24.3.0-0
mamba activate /home/abportillo/.conda/envs/mamba_abner_BC

# get sample_manifest file
rm -f sample_manifest.txt
find . -maxdepth 1 -name "*.bam"  | while read file ; do xbase=$(basename $file); sample_name=${xbase/_sorted_nr_sorted.bam/} ; echo -e "${sample_name}tshortt${xbase}" >> sample_manifest.txt; done
find . -maxdepth 1 -name "*.tab"  | while read file ; do xbase=$(basename $file); sample_name=${xbase/_SJ.out.tab/} ; echo -e "${sample_name}tSJt${xbase}" >> sample_manifest.txt; done
find . -maxdepth 1 -name "*.gtf"  | while read file ; do xbase=$(basename $file); sample_name=${xbase/.gtf/} ; echo -e "${sample_name}tgtft${xbase}" >> sample_manifest.txt; done

# reset all folders if error log
# teprof3 -f sample_manifest.txt --reset


# run teprof3:
# Assemble step (-as) \ # run Process assemble step (-ps) \ # run Filter transcripts (-fs) \\ 
# run Mega assembly (--tacothread) \ # run Quantification (-qs)
teprof3 -f sample_manifest.txt -s 10 -am 1 -as 10 --assemblethread 4 --assemblestrand 1 --assemblejunctionread 2 \
  -ps 10 -pt 0.5 -ptn 100 \
  -fs 10 -fm 1 --filterintronretention 3 \
  --filtermonoexontpm 1 --filterdownstreammate 2 --filterratio 0.5 \
  --tacothread 10 \
  -qs 10 --quansamplenumbercon 100 --quanmode 1 --quanreadlength 75

mamba deactivate  "
  } > "${outdir}/${folder}_runTEProf3.sh"

done < "$folder_list"