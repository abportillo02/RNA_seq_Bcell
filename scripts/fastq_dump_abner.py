
import os
import csv

# Input file with SRR list (1 per line or CSV)
csv_file = "/home/abportillo/github_repo/RNA_seq_Bcell/scripts/SRR_Acc_List.txt"
script_dir = "/home/abportillo/github_repo/RNA_seq_Bcell/scripts/fastq_dump_abner"
outdir = "/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell"
logdir = "/home/abportillo/github_repo/RNA_seq_Bcell/scripts/fastq_dump_abner/logs"

# Make output folders
os.makedirs(script_dir, exist_ok=True)
os.makedirs(outdir, exist_ok=True)
os.makedirs(logdir, exist_ok=True)

with open(csv_file) as f:
    for line in f:
        srr = line.strip()
        if not srr:
            continue

        script_path = os.path.join(script_dir, f"{srr}_download.sh")
        log_file = os.path.join(logdir, f"{srr}.log")

        with open(script_path, "w") as f_out:
            f_out.write(f"""#!/bin/bash
#SBATCH --job-name={srr}_dl
#SBATCH --output={log_file}
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -p all
#SBATCH --mem=50G
#SBATCH --time=03:00:00

module load SRA-Toolkit

echo "Downloading {srr}"
fasterq-dump {srr} -O {outdir} --threads 4
""")

        os.chmod(script_path, 0o755)
        print(f"Created: {script_path}")
