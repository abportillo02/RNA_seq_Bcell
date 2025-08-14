
     
import sys 
import os 

sample =sys.argv[1]
outpath = os.path.join("/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess"
"counts_tx", sample)

outPath = sys.argv[2]
outPath = os.path.join(f"{outPath}", sample)
inputPath = sys.argv[3]

dict = {}
input = open("/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/hg38_p14/Gencode_TE_transcripts.bed", "r")
input = open(f"[inputPath]", "r")

for line in input:
    line = line.strip()
    line = line.split()
    dict[line[3]] = line[4]

input.close() 

#get tx counts intfo for TE transcripts 
input = open(f"{outPath}/{sample}_Gencode_transcripts_ballgown.gtf", "r"))
output = open(f"{outPath}/{sample}_Gencode_TEs.bed", "w")

for line in input:
    line = line.strip()
    line = line.split()
    if line[2] == 'transcript':
        tx_id = line[11].split('"')[1]
        if tx_id in dict.keys():
            output.write(line[0] + '\t' + line[3] + '\t' + line[4] + '\t' +
                         line[9].split('"')[1] + '\t' + tx_id + '\t' +
                         line[6] + "\t" + dict[tx_id] + "\t" +
                         line[-5].split('"')[1] + '\t' + line[-3].split('"')[1] + '\t' +
                         line[-1].split('"')[1] + '\n')
