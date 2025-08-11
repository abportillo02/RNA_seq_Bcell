# get transcripts counts stats, Sep5_2024
# In round4 of counts stats, we want to make this R script generalized for 
# each sample or at list each test folder
# tetx_df_split <- tetx_df %>%
# 1. we first calculate the general numbers: aligned reads, number of txs
# and unmber of TE txst = str_split(`repeat`, pattern = "/", n=6)) %>%
# 2. we get long matrix/df for tx and TE_tx counts numbers, with ids as columns
# and tx names as rowslit, function(x) x[1]),
# 3. we get long matrix/df for tx and TE_tx FPKM numbers, with ids as columns
# and tx names as rowsx) >= 6, paste0(x[2], "/", x[3]), x[2])}),
# 4. we get long matrix/df for tx and TE_tx TPM numbers, with ids as columns
# and tx names as rowsly(split, function(x) {x[(length(x)-1)]}),
# 5. do DE analysis for Tx and TE_tx ) {x[(length(x))]})
# 6. get 5'UTR sequence for high expressed TE_tx
# 7. perform motif enrichment analysis

# laod packages
library(tidyverse)%
library(magrittr)(repFamily == "L1") %>%
library(stringr)e(group = ifelse(repName=="L1HS", "L1HS",
library(rtracklayer)             ifelse(repName=="L1PA2", "L1PA2",
#                                       ifelse(repName=="L1PA3", "L1PA3", "others")))) %>%
# getNum function= sampleID))+
#   geom_bar(aes(fill = group), position = position_stack(reverse = TRUE)) +
getNum <- function(path, pattern){
  # get path and sample names
  counts_path <- list.files(path, recursive = TRUE, pattern = pattern,
                            full.names = TRUE)
  sample_names <- str_remove_all(counts_path, "..+gene_counts/")
  sample_names <- str_remove_all(sample_names, "/..+")L1PA")) %>% 
    dplyr::mutate(group =
                    ifelse(!(repName%in%
  # load data                  c("L1HS", "L1PA2","L1PA3","L1PA4",
  df_list <- list()              "L1PA5", "L1PA6", "L1PA7", "L1PA8", "L1PA8A")),
  for (i in 1:length(sample_names)) {repName)) %>%
    if (str_detect(counts_path[i], "TEs.bed")) {
      # load TEs.bed including FPKM and TPM
      data <- read.table(counts_path[i], = position_stack(reverse = TRUE)) +
                         header = FALSE, sep = "\t",
                         stringsAsFactors = FALSE)
      colnames(data) <- c("seqnames", "start", "end",
                          "gene_id", "transcript_id", "strand", "repeat",
                          "cov", "FPKM", "TPM")
      data$cov = as.numeric(data$cov)
      data$FPKM = as.numeric(data$FPKM)
      data$TPM = as.numeric(data$TPM)
      
    }else if (str_detect(counts_path[i], "filter.gtf")) {
      # load tx.gtf including FPKM and TPM
      data <- import(counts_path[i])
      data <- as.data.frame(data, stringsAsFactors = FALSE)
      data$cov = as.numeric(data$cov)
      data$FPKM = as.numeric(data$FPKM)
      data$TPM = as.numeric(data$TPM)
      
    }else if (str_detect(counts_path[i], "alignedReads.txt")) {
      # load aligned reads number
      data <- read.table(file = counts_path[i], stringsAsFactors = FALSE)
      colnames(data) <- "num_aligned_reads"
      
    }else if (str_detect(counts_path[i], "transcript_count_matrix.csv")) {
      # load csv files
      data <- read_csv(counts_path[i])
      
    }else if (str_detect(counts_path[i], "gene_count_matrix.csv")) {
      data <- read_csv(counts_path[i])
      
    }else if (str_detect(counts_path[i], "TEtxs_count_matrix.csv")) {
      data <- read.table(file = counts_path[i], stringsAsFactors = FALSE)
      colnames(data) <- c("transcript_id", "counts", "repeat")
    }
    
    # add data df to df_list
    df_list[[i]] <- as.data.frame(data)
  }
  
  # name df_list
  names(df_list) <- paste0(sample_names, "_", refGenome)
  
  # retun the df_list
  return(df_list)
}


# # use getNum function to get df lists for 4 refGene test:#######################
# refG <-  c("hg38_p14")


## For TE get df for all number of transcripts##################################
### total num of TE tx and FPKM and TPM
num_tetx_df <- data.frame()
tetx_df <- data.frame()
# for (j in 1:4) {
df_list <-
  getNum(path = paste0("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/Abner/DCIS/project/counts_tx/"),
         pattern = "TEs.bed")

# bind list for TEtx numbers
tetx_count_list <- lapply(df_list, function(x) {nrow(x[x$cov>0,])})
rbind_df_list <- do.call("rbind", tetx_count_list)
rbind_df_list <- as.data.frame(rbind_df_list, stringsAsFactors = FALSE)
colnames(rbind_df_list) <- "num_TEtx"
num_tetx_df <- rbind(num_tetx_df, rbind_df_list)


# bind list for TEtx FPKM and TPM matrix
tetx_list <- lapply(df_list, function(x) {x[x$cov>0,]})
## add id column
tetx_list <- lapply(names(tetx_list), function(x) {
  df <- tetx_list[[x]]
  df$sampleID <- x
  return(df)
})
## rbind for a list
rbind_tetx_list <- do.call("rbind", tetx_list)
## rbind all
tetx_df <- rbind(tetx_df, rbind_tetx_list)
# }



# tetx_df_split <- tetx_df %>%
#   dplyr::mutate(sampleID =str_replace_all(sampleID, "_hg..+", "")) %>%
#   dplyr::mutate(split = str_split(`repeat`, pattern = "/", n=6)) %>%
#   dplyr::mutate(
#     gene = sapply(split, function(x) x[1]),
#     repName = sapply(split, function(x) {
#       ifelse(length(x) >= 6, paste0(x[2], "/", x[3]), x[2])}),
#     repClass = sapply(split, function(x) {x[(length(x)-2)]}),
#     repFamily = sapply(split, function(x) {x[(length(x)-1)]}),
#     teID = sapply(split, function(x) {x[(length(x))]})
#   ) %>%
#   select(-split)

# # plot L1s
# tetx_df_split %>%
#   dplyr::filter(repFamily == "L1") %>%
#   dplyr::mutate(group = ifelse(repName=="L1HS", "L1HS",
#                                ifelse(repName=="L1PA2", "L1PA2",
#                                       ifelse(repName=="L1PA3", "L1PA3", "others")))) %>%
#   ggplot(aes(y = sampleID))+
#   geom_bar(aes(fill = group), position = position_stack(reverse = TRUE)) +
#   theme(legend.position = "top")

# # plot younger L1s
# tetx_df_split %>%
#   dplyr::filter(repFamily == "L1") %>%
#   # dplyr::filter(str_detect(repName, pattern="L1HS|L1PA")) %>% 
#   dplyr::mutate(group =
#                   ifelse(!(repName%in%
#                              c("L1HS", "L1PA2","L1PA3","L1PA4",
#                                "L1PA5", "L1PA6", "L1PA7", "L1PA8", "L1PA8A")),
#                          "others", repName)) %>%
#   dplyr::filter(group!="others") %>%
#   ggplot(aes(y = sampleID))+
#   geom_bar(aes(fill = group), position = position_stack(reverse = TRUE)) +
#   theme(legend.position = "top")
