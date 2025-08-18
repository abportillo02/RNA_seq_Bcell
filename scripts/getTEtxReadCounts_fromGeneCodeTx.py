import sys 
import os 

# /home/abportillo/github_repo/RNA_seq_Bcell/scripts/sampleNames.txt
sample = sys.argv[1]
# outPath = os.path.join("/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess"
                    #    "/counts_tx", sample)
outPath = sys.argv[2]
# outPath = os.path.join(f"{outPath}", f"{sample}")
inputPath =sys.argv[3]

dict = {}
# input = open("/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/hg38_p14/Gencode_TE_transcripts.bed", "r")
input = open(f"{inputPath}", "r")



for line in input:
    line = line.strip()
    line = line.split()
    dict[line[3]] = line[4]

input.close()

# get py script for TE transcripts


# Read sample names from your file
with open(sample, "r") as f:
    samples = [line.strip() for line in f]

# Loop over samples
for s in samples:
    count_matrix1 = os.path.join(outPath, s, f"{s}_transcript_count_matrix.csv")
    count_matrix2 = os.path.join(outPath, s, f"{s}_TEtxs_count_matrix.csv")
    
    # skip if gtf file doesn't exist
    if not os.path.exists(count_matrix1):
        print(f"Warning: {count_matrix1} not found, skipping.")
        continue

    input = open(count_matrix1, "r")
    output = open(count_matrix2, "w")
    
    next(input)
    for line in input:
        line = line.strip()
        tx_id = line.split(",")[0]
        tx_counts = line.split(",")[1]
        if tx_id in dict.keys():
            output.write(tx_id + '\t' + tx_counts + '\t' + dict[tx_id] + '\n')