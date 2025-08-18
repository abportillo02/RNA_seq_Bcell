import sys
import os

def parse_gtf_attributes(attr_str):
    """Parse GTF attribute column into a dictionary."""
    attrs = {}
    for attr in attr_str.strip().split(';'):
        if attr.strip():
            try:
                key, value = attr.strip().split(' ', 1)
                attrs[key] = value.strip('"')
            except ValueError:
                print(f"Malformed attribute: {attr}")
    return attrs

# Input arguments
sample_file = sys.argv[1]      # Path to sampleNames.txt
outPath = sys.argv[2]          # Output directory
inputPath = sys.argv[3]        # BED file with TE annotations

# Load BED entries into a coordinate-based dictionary
coord_dict = {}
with open(inputPath, "r") as input:
    for line in input:
        line = line.strip().split()
        if len(line) >= 5:
            chrom = line[0]
            start = line[1]
            end = line[2]
            key = (chrom, start, end)
            coord_dict[key] = {
                'name': line[3],
                'annotation': line[4]
            }
        else:
            print(f"Skipping malformed BED line: {line}")

print(f"Loaded {len(coord_dict)} coordinate entries from BED file.")

# Read sample names
with open(sample_file, "r") as f:
    samples = [line.strip() for line in f]

print(f"Found {len(samples)} samples to process.")

# Process each sample
for s in samples:
    gtf_file = "/home/abportillo/github_repo/TEProf3/reference/HERV.gtf"
    bed_dir = os.path.join(outPath, s)
    bed_file = os.path.join(bed_dir, f"{s}_Gencode_HERVs.bed")


    os.makedirs(bed_dir, exist_ok=True)


    if not os.path.exists(gtf_file):
        print(f"Warning: {gtf_file} not found, skipping {s}.")
        continue

    print(f"Processing sample: {s}")
    print(f"Writing to: {bed_file}")

    written_lines = 0

    with open(gtf_file, "r") as input, open(bed_file, "w") as output:
        for line in input:
            if line.startswith("#") or line.strip() == "":
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:

                continue

            feature_type = fields[2]
            if feature_type != "exon":
                continue

            chrom = fields[0]
            start = fields[3]
            end = fields[4]
            key = (chrom, start, end)

            attrs = parse_gtf_attributes(fields[8])
            if key in coord_dict:
                entry = coord_dict[key]
                output.write('\t'.join([
                    chrom,
                    start,
                    end,
                    entry['name'],
                    entry['annotation'],
                    fields[6],  # strand from GTF
                    attrs.get("repName", "."),
                    attrs.get("repClass", "."),
                    attrs.get("repFamily", ".")
                ]) + '\n')
                written_lines += 1


    print(f"Finished sample {s}: {written_lines} lines written.")
