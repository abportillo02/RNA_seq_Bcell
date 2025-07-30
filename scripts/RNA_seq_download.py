import sys
import os
import csv

if len(sys.argv) != 3:
    print("Usage: python download_RNAseq_fastq.py <metadata_csv_file> <output_dir>")
    sys.exit(1)

csv_file = sys.argv[1]
base_out_dir = sys.argv[2]

with open(csv_file, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        sampleSRR = row.get("Run")

        if not sampleSRR:
            print(f"Skipping row with missing SRR:{row}")
            continue

        srr_prefix = sampleSRR[:6]
        if len(sampleSRR) <= 9:
            prefix = f"{srr_prefix}"
        else:
            prefix = f"{srr_prefix}/0{sampleSRR[-2:]}"
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
conda activate /home/abportillo/miniconda3/envs/coh

/home/abportillo/miniconda3/envs/coh/bin/ascp -QT -l 300m -P 33001 -k 1 \\
 -i /home/abportillo/miniconda3/envs/coh/etc/aspera/aspera_bypass_rsa.pem \\
 era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/{prefix}/{sampleSRR}/{sampleSRR}_1.fastq.gz \\
 {sample_out_dir}

/home/abportillo/miniconda3/envs/coh/bin/ascp -QT -l 300m -P 33001 -k 1 \\
 -i /home/abportillo/miniconda3/envs/coh/etc/aspera/aspera_bypass_rsa.pem \\
 era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/{prefix}/{sampleSRR}/{sampleSRR}_2.fastq.gz \\
 {sample_out_dir}

conda deactivate
echo "Download completed for {sampleSRR}"
""")

        print(f"Generated script for {sampleSRR}: {script_path}")

        
