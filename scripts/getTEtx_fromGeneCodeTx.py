
import os

# Fixed paths
bed_path = "/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/hg38_p14/Genecode_TE_transcripts.bed"
counts_base = "/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/counts_tx"
output_base = "/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/TE_read_counts"
sample_list_file = "/home/abportillo/github_repo/RNA_seq_Bcell/scripts/sampleNames.txt"  

# Load BED reference into dictionary
bed_dict = {}
with open(bed_path, "r") as bed_file:
    for line in bed_file:
        parts = line.strip().split()
        bed_dict[parts[3]] = parts[4]

# Read sample names
with open(sample_list_file, "r") as f:
    samples = [line.strip() for line in f if line.strip()]

for sample in samples:
    print(f"Processing sample: {sample}")

    gtf_path = os.path.join(counts_base, sample, f"{sample}_Gencode_transcripts_ballgown.gtf")
    out_dir = output_base
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{sample}_Gencode_TEs.bed")

    try:
        with open(gtf_path, "r") as gtf_file, open(out_path, "w") as out_file:
            for line in gtf_file:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split()
                if len(fields) < 12:
                    continue  # skip malformed lines
                
                if fields[2] == "transcript":
                    # Extract transcript ID from attribute field (12th field)
                    # Usually something like: transcript_id "XYZ"; ...
                    attr = fields[11]
                    tx_id = attr.split('"')[1]
                    if tx_id in bed_dict:
                        # Write fields similar to your original script
                        # line[0] = chrom, [3]=start, [4]=end, [9] is attribute, [6]=strand
                        out_file.write(
                            f"{fields[0]}\t{fields[3]}\t{fields[4]}\t"
                            f"{fields[9].split('\"')[1]}\t{tx_id}\t{fields[6]}\t"
                            f"{bed_dict[tx_id]}\t{fields[-5].split('\"')[1]}\t"
                            f"{fields[-3].split('\"')[1]}\t{fields[-1].split('\"')[1]}\n"
                        )
        print(f"Written filtered BED to: {out_path}")
    except FileNotFoundError:
        print(f"Warning: GTF file not found for sample {sample} at {gtf_path}")
