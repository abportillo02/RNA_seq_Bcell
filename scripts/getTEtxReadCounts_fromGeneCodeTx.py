
import os

# Paths fixed for your environment
bed_path = "/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/hg38_p14/Genecode_TE_transcripts.bed"
output_base = "/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/TE_read_counts"
counts_base = "/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/counts_tx"
sample_list_file = "/home/abportillo/github_repo/RNA_seq_Bcell/scripts/sampleNames.txt"  

# Load the BED reference into a dictionary once
bed_dict = {}
with open(bed_path, "r") as bed_file:
    for line in bed_file:
        parts = line.strip().split()
        bed_dict[parts[3]] = parts[4]  # transcriptID -> some value

# Read sample names from file
with open(sample_list_file, "r") as f:
    samples = [line.strip() for line in f if line.strip()]

# Process each sample
for sample in samples:
    print(f"Processing sample: {sample}")
    
    sample_counts_path = os.path.join(counts_base, sample, f"{sample}_transcript_count_matrix.csv")
    sample_out_dir = output_base
    os.makedirs(sample_out_dir, exist_ok=True)
    sample_out_path = os.path.join(sample_out_dir, f"{sample}_TEtxs_count_matrix.csv")
    
    try:
        with open(sample_counts_path, "r") as counts_file, open(sample_out_path, "w") as out_file:
            header = next(counts_file)  # skip header
            
            for line in counts_file:
                line = line.strip()
                if not line:
                    continue
                tx_id, tx_counts = line.split(",")[:2]
                if tx_id in bed_dict:
                    out_file.write(f"{tx_id}\t{tx_counts}\t{bed_dict[tx_id]}\n")
                    
        print(f"Output written to: {sample_out_path}")
        
    except FileNotFoundError:
        print(f"Warning: transcript count matrix file not found for sample {sample}: {sample_counts_path}")
