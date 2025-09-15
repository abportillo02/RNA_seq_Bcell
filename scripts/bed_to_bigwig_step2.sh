# Path to your directory with sorted BED files
dir="/home/abportillo/github_repo/RNA_seq_Bcell/output/chip_exo/KLF4/kzfp_peaks"

# Path to chromosome sizes file
chrom_sizes="hg38.chrom.sizes"  

# Loop through each sorted BED file
for file in "${dir}"/*.sorted.bed; do
  base=$(basename "$file" .sorted.bed)
  bedGraphToBigWig "$file" "$chrom_sizes" "${dir}/${base}.bw"
done
