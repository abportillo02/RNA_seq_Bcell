 


.libPaths()(gene_anno, by = "gene_id")
if (!requireNamespace("BiocManager", quietly = TRUE))_TMM.tsv"))
  install.packages("BiocManager")
lcpm_df <- read_tsv(file = paste0(outDir,"/lcpm_TMM.tsv"))
# use count matrix from feature counts
# get gene counts stats
outDir <- "/xx/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_allo_round2_filter/counts_tx_round4"
# plot
# load packages
library(tidyverse)ne_name %in% c("ZNF141", "MPHOSPH8", 
library(magrittr)                "ATF7IP", "SIRT1",
library(stringr)                 "TASOR", "MORC2", "PPHLN1", "SETDB1", 
library(rtracklayer)             "TASOR2", "IRF2",
library(ggplot2)                 "DNMT1", "DNMT3A")) %>% 
library(ggrepel)
  ggplot(aes(x=sample, y=lcpm))+
# getNum function "identity")+
  labs(title = "HUSH complex genes", x="")+
getNum <- function(path, pattern){"free_y")+
  # get path and sample namesnscript_id, scales = "free_y")+
  counts_path <- list.files(path, recursive = TRUE, pattern = pattern, 
                            full.names = TRUE)
  sample_names <- str_remove_all(basename(counts_path), pattern = "_..+")
  pm_df %>% 
  # get ref. genome versionin% c("ZNF141", "ZNF382", "ZNF429", # young L1 ZNF
  if (str_detect(counts_path[1], "hg19")) { # stem cell ZNF
    refGenome = "hg19"           "NFKB1", "JAK1", "IFIH1", "IRF1", #NFKB(NFKB -> proinflammatory cytokines), IFN(IRF3 -> type I IFNs)
  }else if (str_detect(counts_path[1], "hg38/")) {, #EMT
    refGenome = "hg38_p0"        "E2F1", "E2F2","E2F5", "E2F7", "E2F8", "RB1", "CDK1", #E2F
  }else if (str_detect(counts_path[1], "hg38_seqset/")) {"EED"# PRC2, EZH2 bind to young L1s
    refGenome = "hg38_seqset"    )) %>% 
  }else if (str_detect(counts_path[1], "hg38_p14")) {
    refGenome = "hg38_p14"y')+
  }abs(title = "key genes (ZNF, NFKB, EMT, E2F and PRC)", x="")+
  facet_wrap(~gene_name)+
  # load datatext.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  df_list <- list()
  for (i in 1:length(sample_names)) {
    if (str_detect(counts_path[i], "TEs.bed")) {>%
      # load TEs.bed including FPKM and TPM "CTCF", "HCFC1", "TP53", "THAP11)) %>%
      data <- read.table(counts_path[i], 
                         header = FALSE, sep = "\t", 
                         stringsAsFactors = FALSE)
      colnames(data) <- c("seqnames", "start", "end", 
                          "gene_id", "transcript_id", "strand", "repeat", 
                          "cov", "FPKM", "TPM")
      data$cov = as.numeric(data$cov)LPHA_RESPONSE and HALLMARK_INTERFERON_GAMMA_RESPONSE and KEGG_MEDICUS_REFERENCE_CGAS_STING_SIGNALING_PATHWAY
      data$FPKM = as.numeric(data$FPKM)
      data$TPM = as.numeric(data$TPM)
lcpm_df %>% 
    }else if (str_detect(counts_path[i], "ballgown.gtf")) {
      # load tx.gtf including FPKM and TPM, "FOS", "FOSB", "DUSP1", "SAFB", "SAFB2",
      data <- import(counts_path[i])EM-107", "JUN", "NR4A1", "EGR1", "RGS1",
      data <- as.data.frame(data, stringsAsFactors = FALSE), "ZNF860",
      data <- data[data$type == "transcript",] 
      data$cov = as.numeric(data$cov)
      data$FPKM = as.numeric(data$FPKM)
      data$TPM = as.numeric(data$TPM)x="")+
  facet_wrap(~gene_name, scales = "free_y")+
    }else if (str_detect(counts_path[i], "alignedReads.txt")) {
      # load aligned reads numbert(angle = 90, vjust = 0.5, hjust = 1))
      data <- read.table(file = counts_path[i], stringsAsFactors = FALSE)
      colnames(data) <- "num_aligned_reads"
    }else if (str_detect(counts_path[i], "gene_counts.txt")) {
      # load aligned reads number"ZNF460", "EGR1", "EGR2", "KLF9", "KLF6",
      data <- read.table(file = counts_path[i], stringsAsFactors = FALSE, 
                         header = TRUE)
      colnames(data) <- c("Geneid", "Chr","Start", "End", "Strand",
                          "Length", "crick", "watson")
      data <- data %>% dplyr::mutate(totalCount = crick + watson)
    }else if (str_detect(counts_path[i], "transcript_count_matrix.csv")) {
      # load csv fileselement_text(angle = 90, vjust = 0.5, hjust = 1))
      data <- read_csv(counts_path[i])
      f %>% 
    }else if (str_detect(counts_path[i], "gene_count_matrix.csv")) {
      data <- read_csv(counts_path[i]), 
                                 "ZNF318", "SIX1", "TAZ", "GTF2A2", "GTF2E2",
    }else if (str_detect(counts_path[i], "TEtxs_count_matrix.csv")) {)) %>% 
      data <- read.table(file = counts_path[i], stringsAsFactors = FALSE)
      colnames(data) <- c("transcript_id", "counts", "repeat")
    }s(title = "HUSH complex genes", x="")+
    cet_wrap(~gene_name, scales = "free_y")+
    # add data df to df_listanscript_id, scales = "free_y")+
    df_list[[i]] <- as.data.frame(data)e = 90, vjust = 0.5, hjust = 1))
  }
  
  # name df_listes[order(deGenes$FDR),] %>% 
  names(df_list) <- paste0(sample_names, "_", refGenome)("gene_name"="gene_id"))
  
  # retun the df_list
  return(df_list)
} dplyr::filter(gene_id %in% deGenes$gene_name[1:20]) %>%
  ggplot(aes(x=sample, y=lcpm))+
  geom_bar(stat = 'identity')+
## use getNum FUN to get all transcripts #######################################
  facet_wrap(~gene_name, scales = "free_y")+
df_list <- getNum(path = paste0("/xx/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_allo_round2_filter/counts_tx_round4"),
                  pattern = "gene_counts.txt$")vjust = 0.5, hjust = 1))

deGenes_BMPChigh <- deGenes[deGenes$logFC <= -2,]

# bind list for gene numbers
g_count_list <- lapply(df_list, function(x) {nrow(x[x[,9]>0,])})%>%
num_totalG <- do.call("rbind", g_count_list)
num_totalG <- as.data.frame(num_totalG, stringsAsFactors = FALSE)
colnames(num_totalG) <- "num_G"", x="")+
  facet_wrap(~gene_name, scales = "free_y")+
num_totalG %>% rownames_to_column(var = "samples") %>% free_y")+
  ggplot(aes(x=samples, y=num_G))+(angle = 90, vjust = 0.5, hjust = 1))
  geom_bar(stat = 'identity')+
  labs(title = "detected genes", x="", y = "number of total genes")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
#GO
library(clusterProfiler)
# bind list for gene matrix
gene_list <- lapply(df_list, function(x) {x[x[,9]>0,]})id, "\\.*", "")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
## add id column for each list elementgnc_symbol", "ensembl_gene_id", "entrezgene_id"),
gene_list <- lapply(names(gene_list), function(x) {sion",
    df <- gene_list[[x]]ues = deGenes_BMPChigh$gene_name,
    df$sampleID <- x mart = mart)
    return(df)
    })O <- enrichGO(unique(valid_genes$ensembl_gene_id), OrgDb= 'org.Hs.eg.db', ont="MF", 
#                   keyType = "ENSEMBL", pAdjustMethod = "fdr", pvalueCutoff = 0.25)
## rbind for a list
gene_df <- do.call("rbind", gene_list)
up_GO <- enrichGO(unique(valid_genes$ensembl_gene_id),OrgDb= 'org.Hs.eg.db',ont="BP",  
# load gtf annotation files "ENSEMBL", pAdjustMethod = "fdr", pvalueCutoff = 0.25)
data <- import("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/teAnno_round3/filtered_gencode_v46_chr_patch_hapl_scaff_annotation.gtf")
data <- as.data.frame(data, stringsAsFactors = FALSE)
data <- data[data$type == "transcript",]
# get unique gene_id and gene_name$geneID
gene_anno <- data %>% st(strsplit(preRNA_genes, split = "/"))
  dplyr::select(gene_id, gene_name) %>% unique()reRNA_genes,]

# DE analysis
library(DESeq2)richKEGG(valid_genes$ensembl_gene_id, keyType = "ENSEMBL", 
library(edgeR)          pvalueCutoff = 0.25, organism = "hsa", use_internal_data = TRUE)
count_matrix <- gene_df %>%
  pivot_wider(names_from = sampleID, values_from = totalCount,9)
              id_cols = Geneid) %>% 
  column_to_rownames(var="Geneid") %>% 
  as.matrix()hR)
count_matrix[is.na(count_matrix)] <- 0
colnames(count_matrix) <- str_replace_all(colnames(count_matrix), "_..+", "")
dbs <- listEnrichrDbs()
# plot MDS
lcpm <- cpm(count_matrix, log=TRUE)nes$hgnc_symbol), "KEGG_2021_Human")
library(RColorBrewer)1]], title=names(enriched)[1], showTerms = 3)
col.group <- c(rep("BMPC",3), rep("PMM",8)) %>% factor(levels = c("BMPC", "PMM"))
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
par(mar =c(5.1,4.1,4.1,2.1)) # bottom, left, top, right
plotMDS(lcpm, top = nrow(lcpm), gene.selection = 'common', 
        labels=colnames(count_matrix), col=col.group)
#Avoid Inf when taking log of 0 
# normalize counts with TMM method
dge <- DGEList(counts = count_matrix, 
               group = str_replace_all(colnames(count_matrix), "[:digit:]", ""),
               remove.zeros = TRUE,) -1*score)
dge <- calcNormFactors(dge, method = "TMM")
normalized_lcpm <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE)
# normalized_cpm <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)

plotMDS(normalized_lcpm, top = nrow(normalized_lcpm), gene.selection = 'common', 
        labels=colnames(count_matrix), col=col.group)
m_df = msigdbr(species = "Homo sapiens", category = "H")
hallmark <- m_df[,c("gs_name", "gene_symbol")]
# # Create a DGEList object=hallmark, pvalueCutoff=0.25, pAdjustMethod="BH")
# dge <- DGEList(counts = counts_matrix)
# count_matrix_int <- apply(count_matrix, 2, as.integer)
# rownames(count_matrix_int) <- rownames(count_matrix)

# DE analysis
sample_info <- data.frame(sampleID = colnames(count_matrix))
sample_info$groups <- str_replace_all(sample_info$sampleID, "[:digit:]", "")LAMMATORY_RESPONSE")
gseaplot(em, geneSetID="HALLMARK_HYPOXIA", title="HALLMARK_HYPOXIA")
design <- model.matrix(~0 + groups, data=sample_info)ATION", title="HALLMARK_OXIDATIVE_PHOSPHORYLATION")
con <- makeContrasts(group = groupsPMM-groupsBMPC, levels=design) 
v <- voom(dge, design) %>%LLMARK_DNA_REPAIR", title="HALLMARK_DNA_REPAIR")
  lmFit(design) %>%tID="HALLMARK_P53_PATHWAY", title="HALLMARK_P53_PATHWAY")
  contrasts.fit(con) %>%HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", title="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
  eBayes()m, geneSetID="HALLMARK_PI3K_AKT_MTOR_SIGNALING", title="HALLMARK_PI3K_AKT_MTOR_SIGNALING")
dt <- decideTests(v)ID="HALLMARK_IL6_JAK_STAT3_SIGNALING", title="HALLMARK_IL6_JAK_STAT3_SIGNALING")
summary(dt), geneSetID="HALLMARK_INTERFERON_ALPHA_RESPONSE", title="HALLMARK_INTERFERON_ALPHA_RESPONSE")
gseaplot(em, geneSetID="HALLMARK_INTERFERON_GAMMA_RESPONSE", title="HALLMARK_INTERFERON_GAMMA_RESPONSE")
deGenes <- topTreat(v, coef = 1, number = Inf) %>%title="HALLMARK_MYC_TARGETS_V1")
  tibble::rownames_to_column(var = "gene_id") %>% title="HALLMARK_MYC_TARGETS_V2")
  left_join(data, by = "gene_id") %>%
  distinct(gene_name, .keep_all = TRUE) %>%LING_VIA_NFKB", title="HALLMARK_TNFA_SIGNALING_VIA_NFKB")
  dplyr::filter(adj.P.Val < 0.05 & abs(logFC) > 2) %>%
  dplyr::select(gene_name, seqnames, start, end, strand,="HALLMARK_MITOTIC_SPINDLE")
                transcript_biotype=transcript_type,itle="HALLMARK_G2M_CHECKPOINT")
                logFC, AveExpr, FDR = adj.P.Val)itle="HALLMARK_E2F_TARGETS")


dat <- topTable(v, coef = 1, number=Inf) %>%
  tibble::rownames_to_column(var = "gene_id") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(-t, -B, FDR = `adj.P.Val`) %>%
  mutate(Sig = ifelse(FDR < 0.05 & logFC > 2, "Up", 
                      ifelse(FDR < 0.05 & logFC < -2, "Down", "NS"))) %>% O:BP")
  # add gene anno to top table"gene_symbol")]
  dplyr::left_join(gene_anno, by = "gene_id")gory = "H")
hallmark <- m_df[,c("gs_name", "gene_symbol")]
highlight1 <- dplyr::filter(dat, FDR < 0.05 & abs(logFC) > 2) %>% head(10)H")
highlight2 <- dplyr::filter(dat, FDR < 0.05 & abs(logFC) > 2) %>% arrange(desc(logFC)) %>% head(5)
highlight3 <- dplyr::filter(dat, FDR < 0.05 & abs(logFC) > 2) %>% arrange(logFC) %>% head(5)
highlight <- rbind(highlight1, highlight2, highlight3) %>%
  distinct(gene_name, .keep_all = TRUE)
gseaplot(em, geneSetID="HALLMARK_TNFA_SIGNALING_VIA_NFKB", title="HALLMARK_TNFA_SIGNALING_VIA_NFKB")
dat$label <- ifelse(dat$FDR < 0.01, dat$gene_name, NA)HALLMARK_P53_PATHWAY")
gseaplot(em, geneSetID="HALLMARK_E2F_TARGETS", title="HALLMARK_E2F_TARGETS")
ggplot(dat, aes(x = logFC, y = -log(FDR), colour = as.factor(Sig))) +
  geom_point(size = 2) +ing TF regulators of genes
  scale_color_manual(values=c("cornflowerblue", "grey", "coral1"), name = "Significance") +
  geom_text_repel(aes(label = label), max.overlaps = 10, box.padding = 0.5, size = 3) +,]
  geom_vline(xintercept = 0) +
  theme_bw()
m_df = msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD")
# save deGenes and topTable(dat)symbol")]
readr::write_tsv(deGenes, file = paste0(outDir, "/DE_res.tsv"))od="BH")
readr::write_tsv(dat, file = paste0(outDir,"/topTable_res.tsv"))
# em[str_detect(em@result$core_enrichment, "KLF")]$ID
lcpm_df <- normalized_lcpm %>% as.data.frame() %>% 
  rownames_to_column(var = "gene_id") %>% 
  pivot_longer(cols = c(starts_with("BMPC"), starts_with("PMM")), 
               names_to = "sample", 
               values_to = "lcpm") %>% 
  left_join(gene_anno, by = "gene_id")
readr::write_tsv(lcpm_df, file = paste0(outDir,"/lcpm_TMM.tsv"))

lcpm_df <- read_tsv(file = paste0(outDir,"/lcpm_TMM.tsv"))



# plot
lcpm_df %>% 
  dplyr::filter(gene_name %in% c("ZNF141", "MPHOSPH8", 
                                 "ATF7IP", "SIRT1",
                                 "TASOR", "MORC2", "PPHLN1", "SETDB1", 
                                 "TASOR2", "IRF2",
                                 "DNMT1", "DNMT3A")) %>% 
  distinct() %>%
  ggplot(aes(x=sample, y=lcpm))+
  geom_bar(stat = "identity")+
  labs(title = "HUSH complex genes", x="")+
  facet_wrap(~gene_name, scales = "free_y")+
  # facet_wrap(~gene_name+transcript_id, scales = "free_y")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


lcpm_df %>% 
  dplyr::filter(gene_name %in% c("ZNF141", "ZNF382", "ZNF429", # young L1 ZNF
                                 "ZNF483",  # stem cell ZNF
                                 "NFKB1", "JAK1", "IFIH1", "IRF1", #NFKB(NFKB -> proinflammatory cytokines), IFN(IRF3 -> type I IFNs)
                                 "SNAI2", "TWIST1", #EMT
                                 "E2F1", "E2F2","E2F5", "E2F7", "E2F8", "RB1", "CDK1", #E2F
                                 "SUZ12","EZH2", "PRC1", "EED"# PRC2, EZH2 bind to young L1s
                                 )) %>% 
  ggplot(aes(x=sample, y=lcpm))+
  geom_bar(stat = 'identity')+
  labs(title = "key genes (ZNF, NFKB, EMT, E2F and PRC)", x="")+
  facet_wrap(~gene_name)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

lcpm_df %>% 
  dplyr::filter(str_detect(gene_name, "PTEN")) %>%
  # dplyr::filter(gene_name %in% c("UBE3A", "CTCF", "HCFC1", "TP53", "THAP11)) %>%
  ggplot(aes(x=sample, y=lcpm))+
  geom_bar(stat = 'identity')+
  labs(title = "", x="")+
  facet_wrap(~gene_name)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# need to check HALLMARK_INTERFERON_ALPHA_RESPONSE and HALLMARK_INTERFERON_GAMMA_RESPONSE and KEGG_MEDICUS_REFERENCE_CGAS_STING_SIGNALING_PATHWAY
# add p53

lcpm_df %>% 
  dplyr::filter(gene_name %in% c("CD19", "KLF6", 
                                 "HSP90B1", "FOS", "FOSB", "DUSP1", "SAFB", "SAFB2",
                                 "TMEM-107", "JUN", "NR4A1", "EGR1", "RGS1",
                                 "ATF3", "ZNF165", "ZNF215", "ZNF860",
                                 "ZFP36")) %>% 
  ggplot(aes(x=sample, y=lcpm))+
  geom_bar(stat = "identity")+
  labs(title = "HUSH complex genes", x="")+
  facet_wrap(~gene_name, scales = "free_y")+
  # facet_wrap(~gene_name+transcript_id, scales = "free_y")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


lcpm_df %>% 
  dplyr::filter(gene_name %in% c("ZNF460", "EGR1", "EGR2", "KLF9", "KLF6",
                                 "KLF4", "KLF11", "ZEB1", "PRDM4" )) %>% 
  ggplot(aes(x=sample, y=lcpm))+
  geom_bar(stat = "identity")+
  labs(title = "HUSH complex genes", x="")+
  facet_wrap(~gene_name, scales = "free_y")+
  # facet_wrap(~gene_name+transcript_id, scales = "free_y")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

lcpm_df %>% 
  dplyr::filter(gene_name %in% c("YY1", "ELF1", "FOSB", "JUN",
                                 "YY2", 
                                 "ZNF318", "SIX1", "TAZ", "GTF2A2", "GTF2E2",
                                 "CDC5L", "MAML1", "PSMB5", "ATXN7L3")) %>% 
  ggplot(aes(x=sample, y=lcpm))+
  geom_bar(stat = "identity")+
  labs(title = "HUSH complex genes", x="")+
  facet_wrap(~gene_name, scales = "free_y")+
  # facet_wrap(~gene_name+transcript_id, scales = "free_y")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


deGenes <- deGenes[order(deGenes$FDR),] %>% 
  left_join(unique(data[, c("gene_id", "gene_name")]), c("gene_name"="gene_id"))


lcpm_df %>% 
  dplyr::filter(gene_id %in% deGenes$gene_name[1:20]) %>%
  ggplot(aes(x=sample, y=lcpm))+
  geom_bar(stat = 'identity')+
  labs(title = "top genes", x="")+
  facet_wrap(~gene_name, scales = "free_y")+
  # facet_wrap(~ref_gene_name+transcript_id, scales = "free_y")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

deGenes_BMPChigh <- deGenes[deGenes$logFC <= -2,]

lcpm_df %>% 
  dplyr::filter(gene_id %in% deGenes_BMPChigh$gene_name[31:60]) %>%
  ggplot(aes(x=sample, y=lcpm))+
  geom_bar(stat = 'identity')+
  labs(title = "BMPC high genes", x="")+
  facet_wrap(~gene_name, scales = "free_y")+
  # facet_wrap(~ref_gene_name+transcript_id, scales = "free_y")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))



#GO
library(clusterProfiler)
library(biomaRt)
# ensg_ids <- str_replace_all(res_df_sig_BMPChigh$gene_id, "\\.*", "")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
valid_genes <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id"),
                     filters = "ensembl_gene_id_version",
                     values = deGenes_BMPChigh$gene_name,
                     mart = mart)

# up_GO <- enrichGO(unique(valid_genes$ensembl_gene_id), OrgDb= 'org.Hs.eg.db', ont="MF", 
#                   keyType = "ENSEMBL", pAdjustMethod = "fdr", pvalueCutoff = 0.25)


up_GO <- enrichGO(unique(valid_genes$ensembl_gene_id),OrgDb= 'org.Hs.eg.db',ont="BP",  
                  keyType = "ENSEMBL", pAdjustMethod = "fdr", pvalueCutoff = 0.25)
dotplot(up_GO, orderBy="p.adjust",x="p.adjust",  showCategory=15, font.size=10)
#Plot GO graph
goplot(up_GO)
# preRNA_genes <- up_GO@result[4,]$geneID
# preRNA_genes <- unlist(strsplit(preRNA_genes, split = "/"))
# valid_genes[valid_genes$ensembl_gene_id %in% preRNA_genes,]

#KEGG
# up_KEGG <- enrichKEGG(valid_genes$ensembl_gene_id, keyType = "ENSEMBL", 
#                       pvalueCutoff = 0.25, organism = "hsa", use_internal_data = TRUE)
# 
# dotplot(up_KEGG, orderBy="p.adjust",x="p.adjust", font.size=9)


library(enrichR)

setEnrichrSite("Enrichr") #human
dbs <- listEnrichrDbs()
nrow(dbs)
enriched <- enrichr(unique(valid_genes$hgnc_symbol), "KEGG_2021_Human")
plotEnrich(enriched[[1]], title=names(enriched)[1], showTerms = 3)



# GSEA
## rank result by p values
dum <- dat$FDR
#Avoid Inf when taking log of 0 
dum[which(dum==0)] <- 1e-300
score <- -log10(dum)
#It is important to differentiate up and down regulation
score <- ifelse(dat$logFC>=0, score, -1*score)
names(score) <- dat$gene_name.y
score <- sort(score, decreasing = T)


library(msigdbr)
as.data.frame(msigdbr_collections())
m_df = msigdbr(species = "Homo sapiens", category = "H")
hallmark <- m_df[,c("gs_name", "gene_symbol")]
em <- GSEA(score, TERM2GENE=hallmark, pvalueCutoff=0.25, pAdjustMethod="BH")

#Create a dotplot of top 10 pathways
dotplot(em, orderBy="NES", x="NES")


#Create GSEA Plot for one pathway
gseaplot(em, geneSetID="HALLMARK_INFLAMMATORY_RESPONSE", title="HALLMARK_INFLAMMATORY_RESPONSE")
gseaplot(em, geneSetID="HALLMARK_HYPOXIA", title="HALLMARK_HYPOXIA")
gseaplot(em, geneSetID="HALLMARK_OXIDATIVE_PHOSPHORYLATION", title="HALLMARK_OXIDATIVE_PHOSPHORYLATION")

gseaplot(em, geneSetID="HALLMARK_DNA_REPAIR", title="HALLMARK_DNA_REPAIR")
gseaplot(em, geneSetID="HALLMARK_P53_PATHWAY", title="HALLMARK_P53_PATHWAY")
gseaplot(em, geneSetID="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", title="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
gseaplot(em, geneSetID="HALLMARK_PI3K_AKT_MTOR_SIGNALING", title="HALLMARK_PI3K_AKT_MTOR_SIGNALING")
gseaplot(em, geneSetID="HALLMARK_IL6_JAK_STAT3_SIGNALING", title="HALLMARK_IL6_JAK_STAT3_SIGNALING")
gseaplot(em, geneSetID="HALLMARK_INTERFERON_ALPHA_RESPONSE", title="HALLMARK_INTERFERON_ALPHA_RESPONSE")
gseaplot(em, geneSetID="HALLMARK_INTERFERON_GAMMA_RESPONSE", title="HALLMARK_INTERFERON_GAMMA_RESPONSE")
gseaplot(em, geneSetID="HALLMARK_MYC_TARGETS_V1", title="HALLMARK_MYC_TARGETS_V1")
gseaplot(em, geneSetID="HALLMARK_MYC_TARGETS_V2", title="HALLMARK_MYC_TARGETS_V2")

gseaplot(em, geneSetID="HALLMARK_TNFA_SIGNALING_VIA_NFKB", title="HALLMARK_TNFA_SIGNALING_VIA_NFKB")

gseaplot(em, geneSetID="HALLMARK_MITOTIC_SPINDLE", title="HALLMARK_MITOTIC_SPINDLE")
gseaplot(em, geneSetID="HALLMARK_G2M_CHECKPOINT", title="HALLMARK_G2M_CHECKPOINT")
gseaplot(em, geneSetID="HALLMARK_E2F_TARGETS", title="HALLMARK_E2F_TARGETS")



#We can also use logFC to rank genes
score1 <- res_anno$log2FoldChange
names(score1) <- res_anno$gene_name
score1 <- sort(score1, decreasing = T)
# m_df = msigdbr(species = "Homo sapiens", category = "C5", subcategory="GO:BP")
# GO_BP <- m_df[,c("gs_name", "gene_symbol")]
m_df = msigdbr(species = "Homo sapiens", category = "H")
hallmark <- m_df[,c("gs_name", "gene_symbol")]
em <- GSEA(score1, TERM2GENE=hallmark, pvalueCutoff=0.25, pAdjustMethod="BH")
dotplot(em, orderBy="qvalue", x="NES", font.size=9)

dotplot(em, orderBy="NES", x="NES", font.size=9)

gseaplot(em, geneSetID="HALLMARK_TNFA_SIGNALING_VIA_NFKB", title="HALLMARK_TNFA_SIGNALING_VIA_NFKB")
gseaplot(em, geneSetID="HALLMARK_P53_PATHWAY", title="HALLMARK_P53_PATHWAY")
gseaplot(em, geneSetID="HALLMARK_E2F_TARGETS", title="HALLMARK_E2F_TARGETS")

#use enrichR for predicting TF regulators of genes
# enriched <- enrichr(unique(res_anno_sig_BMPChigh$gene_name), "TRANSFAC_and_JASPAR_PWMs")
# enriched$TRANSFAC_and_JASPAR_PWMs[grep("KLF", enriched$TRANSFAC_and_JASPAR_PWMs$Term),]

#Use msigdb
m_df = msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD")
TFT <- m_df[,c("gs_name", "gene_symbol")]
em <- GSEA(score, TERM2GENE=TFT, pvalueCutoff=0.25, pAdjustMethod="BH")
dotplot(em, orderBy="NES", x="NES")
# em[str_detect(em@result$core_enrichment, "KLF")]$ID

