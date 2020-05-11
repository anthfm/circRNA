library(PharmacoGx)
library(Biobase)
library(calibrate)
library(forestplot)
library(survcomp)
library(DESeq2)
library(dplyr)
library(stringr)

options(stringsAsFactors = F)

matchToIDTable <- function(ids,tbl, column, returnColumn="unique.cellid") {
  sapply(ids, function(x) {
    myx <- grep(paste0("((///)|^)",Hmisc::escapeRegex(x),"((///)|$)"), tbl[,column])
    if(length(myx) > 1){
      stop("Something went wrong in curating ids, we have multiple matches")
    }
    if(length(myx) == 0){return(NA_character_)}
    return(tbl[myx, returnColumn])
  })
}

#read in current cell annotations
cell_all <- read.csv(file = "~/Desktop/Projects/BHKLAB/ncRNA/rnaseq_meta/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
#gcsi metadata
gcsi <- read.csv(file = "~/Desktop/Projects/BHKLAB/ncRNA/rnaseq_meta/gcsi_rnaseq_meta.csv")
gcsi$cellid <- matchToIDTable(ids=gcsi$Cell_line , tbl = cell_all, column = "GNE.cellid", returnColumn = "unique.cellid")
rownames(gcsi) <- gcsi$alias

##########################
#####DESeq2 for CIRI2#####
##########################

#create gene count matrix (# of exonic junction reads >=2 per gene)
load("/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/rnaseq_meta/Gencode.v33lift37.annotation.RData")

summarizeCIRIMatrix <- function(dir_path){
  
  ciri_files <- list.files(path=dir_path,  pattern = "\\.tsv$", full.names = T)
  gene_matrix <- data.frame(matrix(ncol=length(ciri_files), nrow = length(rownames(features_gene))))
  rownames(gene_matrix) <- rownames(features_gene)
  
  for (f in 1:length(ciri_files)) {
    sample <- read.table(file = ciri_files[f], 
                         sep = '\t', 
                         skip = 1,
                         header = FALSE)
    
    colnames(sample) <- c("circRNA_ID","chr","circRNA_start","circRNA_end","junction_reads", "SM_MS_SMS", "non_junction_reads", "junction_reads_ratio", "circRNA_type", "gene_id", "strand","junction_reads_ID")
    sample_name <- gsub("\\..*","", ciri_files[f])
    sample_name <- gsub(".*/","", sample_name)
    sample <- sample[which(sample$circRNA_type=="exon"),]
    sample <- sample[which(sample$junction_reads >= 2),]
    
    gene_reads <- sample %>% 
      group_by(gene_id) %>% 
      summarise(junction_reads = sum(junction_reads))
    
    gene_reads <- as.data.frame(gene_reads)
    gene_reads <- gene_reads[str_count(gene_reads$gene_id, ',')  <= 1,]
    gene_reads$gene_id <- gsub(",", "", gene_reads$gene_id)
    rownames(gene_reads) <- gene_reads$gene_id
    gene_matrix[rownames(gene_reads), f] <- gene_reads$junction_reads
    names(gene_matrix)[f] <- sample_name
    
  }
  gene_matrix <- gene_matrix[rowSums(is.na(gene_matrix)) != ncol(gene_matrix), ]#remove rows that have NO exp (NA) for any sample
  gene_matrix[is.na(gene_matrix)] <- 0
  return(gene_matrix)
} 

gcsi_ciri_matrix <- summarizeCIRIMatrix(dir_path = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRI2/gCSI/result")
#gcsi_ciri_matrix <- gcsi_ciri_matrix + 1 #DeSeq calcualtes geometric mean, which is an issue with counts that have 0. +1 will not affect results, as log is taken anyway
colData <- data.frame(matrix(nrow=ncol(gcsi_ciri_matrix), ncol=2))
colnames(colData) <- c("type", "dataset")
colData$type <- "cell-line"
colData$dataset <- "gCSI"
rownames(colData) <- colnames(gcsi_ciri_matrix)

dds <- DESeqDataSetFromMatrix(countData = gcsi_ciri_matrix, colData = colData, design = ~ 1)
dds2 <- DESeq2::estimateSizeFactors(dds, type='poscounts') #poscounts estimator deals with a gene with some zeros, by calculating a modified geometric mean by taking the n-th root of the product of the non-zero counts
dds_f <- DESeq(dds2)
gcsi_counts_normalized <-counts(dds_f, normalized = TRUE) 



####calculate condordance index####

colnames(gcsi_counts_normalized) <- gsub("gcsi","", colnames(gcsi_counts_normalized))
colnames(gcsi_counts_normalized) <- gcsi$cellid[match(colnames(gcsi_counts_normalized), rownames(gcsi))]

gCSI <- readRDS("/Users/anthmam/Desktop/Projects/BHKLAB/PSets/gCSI.rds")

sensitivity_drug <- summarizeSensitivityProfiles(pSet = gCSI, sensitivity.measure = "aac_recomputed", cell.lines = colnames(gcsi_counts_normalized), fill.missing = F)

ci <- survcomp::concordance.index(sensitivity_drug["Bortezomib",], surv.time = -gcsi_counts_normalized["ENSG00000078808.17_5",colnames(sensitivity_drug)], surv.event = rep(1,ncol(sensitivity_drug)),outx = F, method="noether")



rnaSeq_gCSI <- t(exprs(summarizeMolecularProfiles(gCSI,mDataType = "Kallisto_0.46.1.rnaseq",fill.missing = F)))

# rpkm <- function(counts, lengths) {
#   rate <- counts / lengths 
#   rate / sum(counts) * 1e6
# }
# 
# tpm <- function(counts, lengths) {
#   rate <- counts / lengths
#   rate / sum(rate) * 1e6
# }
