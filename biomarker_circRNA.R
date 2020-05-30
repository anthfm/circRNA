library(PharmacoGx)
library(Biobase)
library(calibrate)
library(forestplot)
library(survcomp)
library(DESeq2)
library(dplyr)
library(stringr)
library(gimme)
#library(plyr)

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

expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  y <- unique(y)
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

#read in current cell annotations
cell_all <- read.csv(file = "~/Desktop/Projects/BHKLAB/ncRNA/rnaseq_meta/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
#gcsi metadata
gcsi <- read.csv(file = "~/Desktop/Projects/BHKLAB/ncRNA/rnaseq_meta/gcsi_rnaseq_meta.csv")
gcsi$cellid <- matchToIDTable(ids=gcsi$Cell_line , tbl = cell_all, column = "GNE.cellid", returnColumn = "unique.cellid")
rownames(gcsi) <- gcsi$alias

##########################
#####CIRI2 Exp Matrix#####
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
gcsi_ciri_matrix_norm <- log2(gcsi_ciri_matrix + 1) #dont do this for deseq2, leave matrix untouched.

##########################
######CIRI2 DEseq2########
##########################
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


##########################
#########CIRI2 CI#########
##########################


colnames(gcsi_ciri_matrix_norm) <- gsub("gcsi","", colnames(gcsi_ciri_matrix_norm))
colnames(gcsi_ciri_matrix_norm) <- gcsi$cellid[match(colnames(gcsi_ciri_matrix_norm), rownames(gcsi))]




gCSI <- readRDS("/Users/anthmam/Desktop/Projects/BHKLAB/PSets/gCSI.rds")
CCLE <- readRDS("/Users/anthmam/Desktop/Projects/BHKLAB/PSets/CCLE.rds")
CTRPv2 <- readRDS("/Users/anthmam/Desktop/Projects/BHKLAB/PSets/CTRPv2.rds")
GDSC2 <- readRDS("/Users/anthmam/Desktop/Projects/BHKLAB/PSets/GDSC2.rds")

combined_mean <- as.data.frame(rowMeans(gcsi_ciri_matrix_norm[, c("SR", "SR")], na.rm = TRUE))
idx <- which(duplicated(colnames(gcsi_ciri_matrix_norm)))
gcsi_ciri_matrix_norm <- gcsi_ciri_matrix_norm[,-idx]
v <- combined_mean$`rowMeans(gcsi_ciri_matrix_norm[, c("SR", "SR")], na.rm = TRUE)`
gcsi_ciri_matrix_norm$SR <- v
  


sensitivity_drug <- summarizeSensitivityProfiles(pSet = gCSI, sensitivity.measure = "aac_recomputed", fill.missing = F)
#missingSamples <- colnames(gcsi_ciri_matrix_norm)[which(!colnames(gcsi_ciri_matrix_norm) %in% colnames(sensitivity_drug))]

sensitivity_drug_ctrpv2 <- as.data.frame(summarizeSensitivityProfiles(pSet = CTRPv2, sensitivity.measure = "aac_recomputed", cell.lines = colnames(gcsi_ciri_matrix_norm), fill.missing = F))
sensitivity_drug_ccle <- as.data.frame(summarizeSensitivityProfiles(pSet = CCLE, sensitivity.measure = "aac_recomputed", cell.lines = colnames(gcsi_ciri_matrix_norm), fill.missing = F))
sensitivity_drug_gdsc <- as.data.frame(summarizeSensitivityProfiles(pSet = GDSC2, sensitivity.measure = "aac_recomputed", cell.lines = colnames(gcsi_ciri_matrix_norm), fill.missing = F))

#aa <- rbind.fill(sensitivity_drug, sensitivity_drug_ctrpv2, sensitivity_drug_gdsc,sensitivity_drug_ccle)
#aa <- rbind(as.matrix(sensitivity_drug), as.matrix(sensitivity_drug_ctrpv2), as.matrix(sensitivity_drug_gdsc),as.matrix(sensitivity_drug_ccle))

#sensitivity_drug$row.names <- rownames(sensitivity_drug)
#sensitivity_drug_ctrpv2$row.names <- rownames(sensitivity_drug_ctrpv2)
#sensitivity_drug_gdsc$row.names <- rownames(sensitivity_drug_gdsc)
#sensitivity_drug_ccle$row.names <- rownames(sensitivity_drug_ccle)


#ci <- survcomp::concordance.index(sensitivity_drug["Bortezomib",commonSamples], surv.time = unlist(-gcsi_ciri_matrix_norm["ENSG00000078808.17_5", commonSamples]), surv.event = rep(1,length(sensitivity_drug[commonSamples])),outx = F, method="noether")

commonSamples <- intersect(colnames(sensitivity_drug),colnames(gcsi_ciri_matrix_norm))

computeCI <- function(circRNA_normalized_data, sensitivity_data, samples){
  
  genes <- rownames(circRNA_normalized_data)
  drugs <- unique(rownames(sensitivity_data))
  commonSamples <- samples
  combinations <- as.data.frame(expand.grid.unique(genes, drugs, include.equals = TRUE))
  combinations$ci <- NA
  combinations$pvalue <- NA
  colnames(combinations) <- c("gene","drug","ci","pvalue")
  
  for (i in 1:nrow(combinations)){
    print(paste0(i, " out of ", nrow(combinations), "complete"))
    tt <- sensitivity_drug[combinations[,2][i],commonSamples]
    tt[which(is.na(tt))] <- 0 #some sensitivities are NA due to filterNoisyCurve function, which causes error when running CI with survcomp
    ci <- survcomp::concordance.index(tt, surv.time = unlist(-gcsi_ciri_matrix_norm[combinations[,1][i], commonSamples]), surv.event = rep(1,length(sensitivity_drug[commonSamples])),outx = F, method="noether")
    combinations$pvalue[i] <- ci$p.value
    combinations$ci[i] <- ci$c.index
  }
  
  return(combinations)
}

gCSI_sens <- computeCI(circRNA_normalized_data=gcsi_ciri_matrix_norm, sensitivity_data=sensitivity_drug, samples = commonSamples)

gCSI_sens_filtered <- gCSI_sens[which(gCSI_sens$pvalue < 0.05 & gCSI_sens$ci > 0.5),]















ci <- survcomp::concordance.index(tt, surv.time = unlist(-gcsi_ciri_matrix_norm["ENSG00000237973.1_4", commonSamples]), surv.event = rep(1,length(sensitivity_drug[commonSamples])),outx = F, method="noether")

commonSamples <- intersect(colnames(sensitivity_drug),colnames(gcsi_ciri_matrix_norm))
genes <- rownames(gcsi_ciri_matrix_norm)
drugs <- unique(rownames(sensitivity_drug))
combinations <- as.data.frame(expand.grid.unique(genes, drugs, include.equals = TRUE))
combinations$ci <- NA
combinations$pvalue <- NA
colnames(combinations) <- c("gene","drug","ci","pvalue")
for (i in 1:nrow(combinations)){
  tt <- sensitivity_drug[combinations[,2][i],commonSamples]
  tt[which(is.na(tt))] <- 0 #some are NA due to filterNoisyCurve function, which causes error when running CI with survcomp
  ci <- survcomp::concordance.index(tt, surv.time = unlist(-gcsi_ciri_matrix_norm[combinations[,1][i], commonSamples]), surv.event = rep(1,length(sensitivity_drug[commonSamples])),outx = F, method="noether")
  combinations$pvalue[i] <- ci$p.value
  combinations$ci[i] <- ci$c.index
}

combinations_fil <- combinations[which(combinations$pvalue < 0.05 & combinations$ci > 0.5),]




rnaSeq_gCSI <- t(exprs(summarizeMolecularProfiles(gCSI,mDataType = "Kallisto_0.46.1.rnaseq", fill.missing = F)))



# rpkm <- function(counts, lengths) {
#   rate <- counts / lengths 
#   rate / sum(counts) * 1e6
# }
# 
# tpm <- function(counts, lengths) {
#   rate <- counts / lengths
#   rate / sum(rate) * 1e6
# }







summarizeCIRCMatrix <- function(dir_path){
  
  circ_files <- list.files(path=dir_path,  pattern = ".txt", full.names = T, recursive = TRUE)
  gene_matrix <- data.frame(matrix(ncol=length(circ_files), nrow = length(rownames(features_gene))))
  rownames(gene_matrix) <- rownames(features_gene)
  
  for (f in 1:length(circ_files)) {
    sample <- read.table(file = circ_files[f], 
                         sep = '\t', 
                         header = FALSE)
    
    colnames(sample) <- c("chrom","start","end", "name","score","strand","thickStart","thickEnd","itemRgb","exonCount","exonSizes", "exonOffsets", "readNumber", "circType", "geneName", "isoformName", "index", "flankIntron")
    sample_name <- sub(".*/ *(.*?) *.txt.*", "\\1", circ_files[f])
    sample_name <- gsub("_circularRNA_known", "", sample_name)
    sample <- sample[which(sample$circType=="circRNA"),]
    sample <- sample[which(sample$readNumber >= 2),]
    
    #get gene-id from gencode v33lift37 feature_transcript data-frame, as circexplorer2 does not provide it in result files
    
    sample$gene_id <- features_transcript$gene_id[match(sample$isoformName, features_transcript$transcript_id)]
    
    #group sum of read number per gene_id
    gene_reads <- sample %>% 
      group_by(gene_id) %>% 
      summarise(readNumber = sum(readNumber))
    
    gene_reads <- as.data.frame(gene_reads)
    
    rownames(gene_reads) <- gene_reads$gene_id
    gene_matrix[rownames(gene_reads), f] <- gene_reads$readNumber
    names(gene_matrix)[f] <- sample_name
    
  }
  gene_matrix <- gene_matrix[rowSums(is.na(gene_matrix)) != ncol(gene_matrix), ]#remove rows that have NO exp (NA) for any sample
  gene_matrix[is.na(gene_matrix)] <- 0
  return(gene_matrix)
} 


gcsi_test <- summarizeCIRCMatrix(dir_path = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRCexplorer2/gCSI/unmapped_method/annotate")
