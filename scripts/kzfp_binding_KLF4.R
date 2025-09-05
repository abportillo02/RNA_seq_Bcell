# ====== Load Required Libraries ======
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(plyranges)
library(AnnotationHub)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

# ====== Set Working Directory ======
dir <- "/home/abportillo/github_repo/RNA_seq_Bcell/scripts/chip_exo/KLF4"

# ====== Load hg19 to hg38 LiftOver Chain ======
ah.chain <- subset(AnnotationHub(ask = FALSE), 
                   rdataclass == "ChainFile" & species == "Homo sapiens")
query(ah.chain, c("hg19", "hg38"))
# get chain
chain <- ah.chain[ah.chain$title == "hg19ToHg38.over.chain.gz"]
chain <- chain[[1]]

# ====== Get KLF4 Genomic Coordinates (hg38) ======
klf4_entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = "KLF4", keytype = "SYMBOL", columns = "ENTREZID")$ENTREZID
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
klf4_gr <- genes(txdb)[klf4_entrez]

# ====== List ChIP-exo Files for KZFPs ======
exo_files <- list.files(
  path = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/chipExo/GSE78099", 
  pattern = "*_exo.bed.gz", full.names = TRUE
)

# ====== Function to Find KZFPs Binding to KLF4 ======
get_KZFP_binding_to_KLF4 <- function(klf4_gr) {
  all_overlaps <- list()
  
  for (exo_file in exo_files) {
    KZFP_id <- str_remove_all(basename(exo_file), "_peaks.+")
    KZFP_name <- str_remove_all(KZFP_id, ".+_")
    KZFP_name <- str_remove_all(KZFP_name, "-rep.")
    # Load current exo.bed.gz files 
    exo_ranges <- rtracklayer::import(exo_file)
    exo_ranges$KZFP_id <- rep(KZFP_id, length(exo_ranges))
    exo_ranges$KZFP_name <- rep(KZFP_name, length(exo_ranges))

    # Perform liftover to hg38
    exo_ranges_list <- liftOver(exo_ranges, chain)
    # Unlist and check the results
    exo_ranges_hg38 <- unlist(exo_ranges_list)

    # Find overlaps with KLF4
    overlaps <- find_overlaps(klf4_gr, exo_ranges_hg38)
    overlaps_df <- as.data.frame(overlaps)

    if (nrow(overlaps_df) > 0) {
      overlaps_df$KZFP_id <- KZFP_id
      overlaps_df$KZFP_name <- KZFP_name
      all_overlaps[[basename(exo_file)]] <- overlaps_df
    }
  }

  final_overlaps <- bind_rows(all_overlaps, .id = "exo_basename")
  return(final_overlaps)
}

# ====== Run the Function ======
all_KZFPs_binding_KLF4 <- get_KZFP_binding_to_KLF4(klf4_gr)

# ====== Save Results ======
write_tsv(all_KZFPs_binding_KLF4, paste0(dir, "/KZFP_binding_to_KLF4.tsv"))
