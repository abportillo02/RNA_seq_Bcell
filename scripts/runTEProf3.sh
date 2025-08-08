#!/bin/bash

folder=$1
# Define the directory:
while IFS= read -r folder; do
  echo "Creating bash script for folder: $folder"

  # teprof3Dir="/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/Abner/DCIS/project/Transcript_finder/TEProf3_TEdetect"

  # cd ${teprof3Dir}

  outdir=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/Abner/DCIS/project/Transcript_finder/${folder}

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

# Softlink data to 4 folders 

# ls /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/Abner/DCIS/project/*/*_sorted_nr_sorted.bam.bai | sed -n '1,10p' | xargs -I {} ln -s "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/Abner/DCIS/project/*/*_sorted_nr_sorted.bam.bai" /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/Abner/DCIS/project/Transcript_finder/TEProf3_TEdetect1

source /home/abportillo/.bashrc
module load Mamba/24.3.0-0
mamba activate /home/abportillo/.conda/envs/mamba_abner_BC

 
# get sample_manifest file
# rm sample_manifest.txt
# find . -maxdepth 1 -name "*.bam"  | while read file ; do xbase=$(basename $file); sample_name=${xbase/_sorted_nr_sorted.bam/} ; echo -e "${sample_name}\tshort\t${xbase}" >> sample_manifest.txt; done ;
# find . -maxdepth 1 -name "*.tab"  | while read file ; do xbase=$(basename $file); sample_name=${xbase/_SJ.out.tab/} ; echo -e "${sample_name}\tSJ\t${xbase}" >> sample_manifest.txt; done ;
# find . -maxdepth 1 -name "*.gtf"  | while read file ; do xbase=$(basename $file); sample_name=${xbase/.gtf/} ; echo -e "${sample_name}\tgtf\t${xbase}" >> sample_manifest.txt; done

# reset all folders if error log
teprof3 -f sample_manifest.txt --reset


# run teprof3:ZZ
# Assemble step (-as) \ # run Process assemble step (-ps) \ # run Filter transcripts (-fs) \\ 
# run Mega assembly (--tacothread) \ # run Quantification (-qs)
teprof3 -f sample_manifest.txt -s 10 -am 1 -as 10 --assemblethread 4 --assemblestrand 1 --assemblejunctionread 2 \
  -ps 10 -pt 0.5 -ptn 100 \
  -fs 10 -fm 1 --filterintronretention 3 \
  --filtermonoexontpm 1 --filterdownstreammate 2 --filterratio 0.5 \
  --tacothread 10 \
  -qs 10 --quansamplenumbercon 100 --quanmode 1 --quanreadlength 75

mamba deactivate  "
 
  } > ${outdir}/${folder}_runTEProf3.sh
#  cd ${outdir}
#  sbatch ${folder}_runTEProf3.sh
done < ${folder}

# # change .bashrc
# vim ~/.bashrc 
# # add these lines:
# PATH=$PATH:/home/abportillo/.conda/envs/mamba_abner_BC/bin/TEProf3/bin/
# export PATH="/home/abportillo/.conda/envs/taco_env/bin:$PATH"