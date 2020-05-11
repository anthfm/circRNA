library(PharmacoGx)
library(Biobase)
library(calibrate)
library(forestplot)
library(survcomp)
library(DESeq2)
library(dplyr)
library(stringr)

options(stringsAsFactors = F)


sample <- read.table(file = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRI2/gCSI/result/gcsi636438.tsv", 
                     sep = '\t', 
                     skip = 1,
                     header = FALSE)
colnames(sample) <- c("circRNA_ID","chr","circRNA_start","circRNA_end","junction_reads", "SM_MS_SMS", "non_junction_reads", "junction_reads_ratio", "circRNA_type", "gene_id", "strand","junction_reads_ID")
sample <- sample[which(sample$circRNA_type=="exon"),]
sample <- sample[which(sample$junction_reads >= 2),]

sample <- sample[,c(1,2,3,4,5,9,10)]

library(dplyr)


##########################
#####DESeq2 for CIRI2#####
##########################

#create gene count matrix (# of exonic junction reads >=2 per gene)
load("/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/rnaseq_meta/Gencode.v33lift37.annotation.RData")
gene_matrix <- data.frame(matrix(ncol=49, nrow = length(rownames(features_gene))))
rownames(gene_matrix) <- rownames(features_gene)

gene_reads <- sample %>% 
  group_by(gene_id) %>% 
  summarise(junction_reads = sum(junction_reads))

gene_reads <- as.data.frame(gene_reads)
gene_reads <- gene_reads[str_count(gene_reads$gene_id, ',')  <= 1,]
gene_reads$gene_id <- gsub(",", "", gene_reads$gene_id)
rownames(gene_reads) <- gene_reads$gene_id
gene_matrix[rownames(gene_reads), 1] <- c(sample_name, gene_reads$junction_reads)


summarizeCIRIMatrix <- function(dir_path){
  
  ciri_files <- list.files(path=dir_path,  pattern = "\\.tsv$", full.names = T)
  circ_counts_df <- data.frame(matrix(ncol=2, nrow = length(ciri_files)))
  colnames(circ_counts_df) <- c("sample", "count")
  for (f in 1:length(ciri_files)) {
    sample <- read.table(file = ciri_files[f], 
                         sep = '\t', 
                         skip = 1,
                         header = FALSE)
    
    colnames(sample) <- c("circRNA_ID","chr","circRNA_start","circRNA_end","junction_reads", "SM_MS_SMS", "non_junction_reads", "junction_reads_ratio", "circRNA_type", "gene_id", "strand","junction_reads_ID")
    sample_name <- gsub("\\..*","", ciri_files[f])
    sample_name <- gsub(".*/","", sample_name)
    #sample_name <- gsub("gcsi", "", sample_name)
    #sample_name <- gsub("gdsc", "", sample_name)
    
    sample <- sample[which(sample$circRNA_type=="exon"),]
    sample <- sample[which(sample$junction_reads >= 2),]
    count <- sum(sample$junction_reads)
    circ_counts_df[f,] <- c(sample_name,count)
  }
  #rownames(circ_counts_df) <- circ_counts_df$sample
  return(circ_counts_df)
} 









# rpkm <- function(counts, lengths) {
#   rate <- counts / lengths 
#   rate / sum(counts) * 1e6
# }
# 
# tpm <- function(counts, lengths) {
#   rate <- counts / lengths
#   rate / sum(rate) * 1e6
# }
