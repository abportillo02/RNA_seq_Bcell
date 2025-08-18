import sys
import os

def parse_gtf_attributes(attr_str):
    """Parse GTF attribute column into a dictionary."""
    attrs = {}
    for attr in attr_str.strip().split(';'):
        if attr.strip():
            key, value = attr.strip().split(' ', 1)
            attrs[key] = value.strip('"')
    return attrs

# --- Input arguments ---
sample_file = sys.argv[1]      # Path to sampleNames.txt
outPath = sys.argv[2]          # Output directory
inputPath = sys.argv[3]        # BED file with transcript annotations

# --- Load transcript annotations into dictionary ---
dict = {}
with open(inputPath, "r") as input:
    for line in input:
        line = line.strip().split()
        if len(line) > 4:
            dict[line[3]] = line[4]
        else:
            print(f"Skipping malformed line: {line}")

# --- Read sample names ---
with open(sample_file, "r") as f:
    samples = [line.strip() for line in f]

# --- Process each sample ---
for s in samples:
    gtf_file = "/home/abportillo/github_repo/TEProf3/reference/HERV.gtf"
    bed_dir = os.path.join(outPath, s)
    bed_file = os.path.join(bed_dir, f"{s}_Gencode_HERVs.bed")

    # Ensure output directory exists
    os.makedirs(bed_dir, exist_ok=True)

    # Skip if GTF file doesn't exist
    if not os.path.exists(gtf_file):
        print(f"Warning: {gtf_file} not found, skipping.")
        continue

    print(f"Processing sample: {s}")
    print(f"Writing to: {bed_file}")

    with open(gtf_file, "r") as input, open(bed_file, "w") as output:
        for line in input:
            if line.startswith("#") or line.strip() == "":
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type != "transcript":
                continue

            attrs = parse_gtf_attributes(fields[8])
            tx_id = attrs.get("transcript_id")

            if tx_id and tx_id in dict:
                output.write('\t'.join([
                    fields[0],  # chrom
                    fields[3],  # start
                    fields[4],  # end
                    attrs.get("gene_id", "."),  # gene name
                    tx_id,
                    fields[6],  # strand
                    dict[tx_id],
                    attrs.get("repName", "."),
                    attrs.get("repClass", "."),
                    attrs.get("repFamily", ".")
                ]) + '\n')
            else:
                print(f"Transcript ID {tx_id} not found in dict.")
