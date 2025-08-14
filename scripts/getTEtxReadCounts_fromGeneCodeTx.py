import sys 
import os 

sample = sys.argv[1]
outPATH = os.path.join("/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess"
                       "/counts_tx", sample)
outPath = sys.argv[2]
outPath = os.path.join(f"{outPath}", f"{sample}")
input = open(f"{inputPath}", "r")

for line in input:
    line = line.strip()
    line = line.split()
    dict[line[3]] = line[4]

input.close()

# get py script for TE transcripts
input = open(f"{outPath}/{sample}_transcript_count_matrix.csv", "r")
output = open(f"{outPath}/{sample}_TEtxs_count_matrix.csv", "w")

next(input)
for line in input:
    line = line.strip()
    tx_id = line.split(",")[0]
    tx_counts = line.split(",")[1]
    if tx_id in dict.keys():
        output.write(tx_id + '\t' + tx_counts + '\t' + dict[tx_id] + '\n')