import sys
import os
import csv

def sra_prefix(srr):
    """
    Returns the correct SRA FTP/Aspera prefix for a given SRR accession.
    """
    if len(srr) >= 7:
        # e.g. SRR123456 -> SRR123/SRR123456/SRR123456
        return f"{srr[:6]}/{srr}/{srr}"
    else:
        # fallback for shorter SRRs
        return f"{srr[:6]}/{srr}"
        
if len(sys.argv) != 3:
    print("Usage: python download_RNAseq_fastq.py <metadata_csv_file> <output_dir>")
    sys.exit(1)

csv_file = sys.argv[1]
base_out_dir = sys.argv[2]

with open(csv_file) as f:
    reader = csv.DictReader(f)
    for row in reader:
        sampleSRR = row.get("Run")

        if not sampleSRR:
            print(f"Skipping row with missing SRR:{row}")
            continue

        prefix = sra_prefix(sampleSRR)
        print(f"Processing sample: {sampleSRR} with prefix:{prefix}")

        sample_out_dir = os.path.join(base_out_dir, sampleSRR)
        os.makedirs(sample_out_dir, exist_ok=True)

        script_path = os.path.join(sample_out_dir, f"{sampleSRR}_rawData_rna.sh")

        with open(script_path, "w") as script_file:
            script_file.write(f"""#!/bin/bash
#SBATCH --job-name=dl_{sampleSRR}
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -p all
#SBATCH --mem=50G
#SBATCH --time=05:30:00
#SBATCH --output=dl_{sampleSRR}_%j.log

source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC

ascp -QT -l 300m -P 33001 -k 1 \\
 -i /home/abportillo/.ssh/id_rsa.pub.pub\\
 era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/{prefix}_1.fastq.gz \\
 .

ascp -QT -l 300m -P 33001 -k 1 \\
-i /home/abportillo/.ssh/id_rsa.pub.pub\\
 era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/{prefix}_2.fastq.gz \\
.

conda deactivate
echo "Download completed for {sampleSRR}"
""")

        print(f"Generated script for {sampleSRR}: {script_path}")