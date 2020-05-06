##script to read in processed circRNA data (CCLE, GDSC, gCSI, Hansen) from CIRCexplorer2 and CIRI2 for comparison

options(stringsAsFactors = F)

####################################################
####### Determine replicates across datasets #######
####################################################

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

#gcsi
gcsi <- read.csv(file = "~/Desktop/Projects/BHKLAB/ncRNA/rnaseq_meta/gcsi_rnaseq_meta.csv")
gcsi$cellid <- matchToIDTable(ids=gcsi$Cell_line , tbl = cell_all, column = "GNE.cellid", returnColumn = "unique.cellid")
rownames(gcsi) <- gcsi$alias

#ccle
ccle <- read.csv(file = "~/Desktop/Projects/BHKLAB/ncRNA/rnaseq_meta/ccle_rnaseq_meta.csv")
ccle$cellid <- matchToIDTable(ids=ccle$Cell_Line , tbl = cell_all, column = "CCLE.cellid", returnColumn = "unique.cellid")
rownames(ccle) <- ccle$Run

#gdsc
gdsc <- read.csv(file = "~/Desktop/Projects/BHKLAB/ncRNA/rnaseq_meta/gdsc_rnaseq_meta.txt", sep = "\t")
gdsc <- gdsc[which(!gdsc$Comment.SUBMITTED_FILE_NAME. == "15552_5.cram"),]
gdsc$cellid <- matchToIDTable(ids=gdsc$Source.Name, tbl = cell_all, column = "GDSC_rnaseq.cellid", returnColumn = "unique.cellid")
gdsc$files <- gsub(".cram","",gdsc$Comment.SUBMITTED_FILE_NAME.)
rownames(gdsc) <- gdsc$files

#filter metadata -> keep intersected cell lines only
intersected_rnacells <- Reduce(intersect, list(gdsc$cellid, ccle$cellid, gcsi$cellid))

gcsi <- gcsi[which(gcsi$cellid %in% intersected_rnacells),]
ccle <- ccle[which(ccle$cellid %in% intersected_rnacells),]
gdsc <- gdsc[which(gdsc$cellid %in% intersected_rnacells),]


####################################################
######## Summarize # of circRNA's for CIRI2 ########
####################################################

####### read in CIRI2 abundance of circRNA's #######

#reads in CIRI2 abundance data (# of junction reads for each sample)
summarizeCIRI <- function(dir_path){
  
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

gcsi_ciri_counts <- summarizeCIRI(dir_path = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRI2/gCSI/result")
ccle_ciri_counts <- summarizeCIRI(dir_path = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRI2/CCLE/result")
gdsc_ciri_counts <- summarizeCIRI(dir_path = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRI2/GDSC/result")
hansen_ciri_counts <- summarizeCIRI(dir_path = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRI2/hansen/result")
hansen_ciri_matched <- summarizeCIRI(dir_path = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRI2/hansen_match/result")

gcsi_ciri_counts$sample <- gsub("gcsi", "", gcsi_ciri_counts$sample)
gdsc_ciri_counts$sample <- gsub("gdsc", "", gdsc_ciri_counts$sample)

rownames(gcsi_ciri_counts) <- gcsi_ciri_counts$sample
rownames(ccle_ciri_counts) <- ccle_ciri_counts$sample
rownames(gdsc_ciri_counts) <- gdsc_ciri_counts$sample

#average SR technical replicate circRNA count in gCSI
gcsi_ciri_counts$cellid <- gcsi$cellid[match(rownames(gcsi_ciri_counts), rownames(gcsi))]
gcsi_SR_mean <- mean(as.numeric(gcsi_ciri_counts$count[which(gcsi_ciri_counts$cellid=="SR")]))
gcsi_ciri_counts[which(gcsi_ciri_counts$cellid=="SR"),][2,"count"] <- gcsi_SR_mean
gcsi_ciri_counts <- gcsi_ciri_counts[-which(gcsi_ciri_counts$sample=="587641"),]

#order ccle & gdsc cellid by gcsi cellid
ccle_ciri_counts$cellid <- ccle$cellid[match(rownames(ccle_ciri_counts), rownames(ccle))]
ccle_ciri_counts <- ccle_ciri_counts[match(gcsi_ciri_counts$cellid, ccle_ciri_counts$cellid),]

gdsc_ciri_counts$cellid <- gdsc$cellid[match(rownames(gdsc_ciri_counts), rownames(gdsc))]
gdsc_ciri_counts <- gdsc_ciri_counts[match(gcsi_ciri_counts$cellid, gdsc_ciri_counts$cellid),]

ciri_combined <- data.frame("gcsi"=as.numeric(gcsi_ciri_counts$count), "ccle"=as.numeric(ccle_ciri_counts$count), "gdsc"=as.numeric(gdsc_ciri_counts$count))
rownames(ciri_combined) <- gcsi_ciri_counts$cellid
colnames(ciri_combined) <- c("gCSI","CCLE","GDSC")


#violin plot of ccle, gdsc, gcsi circRNA counts
library(ggplot2)
library(reshape2)
df.m <- reshape2::melt(ciri_combined, id.vars = NULL)
pdf(file="/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/circRNA_plots/CIRI2_violin.pdf", width = 10)
ggplot(df.m, aes(x = variable, y = value, fill=variable)) + geom_violin() + 
  stat_summary(
  fun.data = "mean_sdl",  fun.args = list(mult = 1), 
  geom = "pointrange", color = "black")  + scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))  + labs(title="circRNA abundance per dataset (CIRI2)", 
                                                                                                                subtitle="",
                                                                                                                x="Dataset",
                                                                                                                y="circRNA Number") + theme(legend.position="none", plot.title = element_text(hjust = 0.5))
dev.off()


#violin plot of ccle, gdsc, gcsi circRNA counts (log2)
library(ggplot2)
library(reshape2)
df.m <- reshape2::melt(ciri_combined, id.vars = NULL)
#df.m$value <- log2(df.m$value + 1) # Get log2 counts per million - because the data is count data (gene counts) and will contain many 0s we need to add a count of 1 to every value in order to prevent attempting log2(0) from creating errors. This data has no 0's but I will use this to keep consistent when doing gene expression analysis next.
df.m$value <- log2(df.m$value)
pdf(file="/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/circRNA_plots/CIRI2_violin_log2.pdf", width = 10)
ggplot(df.m, aes(x = variable, y = value, fill=variable)) + geom_violin() + 
  stat_summary(
    fun.data = "mean_sdl",  fun.args = list(mult = 1), 
    geom = "pointrange", color = "black")  + scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))  + labs(title="circRNA abundance per dataset (CIRI2)", 
                                                                                                                    subtitle="",
                                                                                                                    x="Dataset",
                                                                                                                    y="Log2(abundance)") + theme(legend.position="none", plot.title = element_text(hjust = 0.5))
dev.off()

#heatmap of ccle, gdsc, gcsi circRNA counts (log2)
library(gplots)
ciri_combined_t <- t(ciri_combined)
ciri_combined_t <- log2(ciri_combined_t)

pdf(file="/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/circRNA_plots/CIRI2_heatmap.pdf", width = 10)
heatmap.2(ciri_combined_t, col=bluered,trace="none", main="CIRI2 circRNA counts (log2)",scale="row", density.info="none", margins=c(12,8), keysize=1)
dev.off()


#boxplot of ccle, gdsc, gcsi circRNA counts (log2)
my_colors <- c( 
  rgb(143,199,74,maxColorValue = 255),
  rgb(242,104,34,maxColorValue = 255), 
  rgb(111,145,202,maxColorValue = 255)
)

pdf(file="/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/circRNA_plots/CIRI2_boxplot_log2.pdf", width = 10)
boxplot(log2(ciri_combined), 
        xlab="Dataset", 
        ylab="Log2(Abundance)",
        las=2,
        col=my_colors)
dev.off()



#filter circRNA counts by at least 1.5-fold enrichment in matching RNAse-R sample
filterCIRI <- function(nonRNAseR_dir, RNAseR_dir){
  #"*neg.tsv$" <- for ribominus
  #"\\.tsv$"<- for poly-A selected data
  ciri_files <- list.files(path=nonRNAseR_dir,  pattern = "*neg.tsv$", full.names = T)
  circ_counts_df <- data.frame(matrix(ncol=2, nrow = length(ciri_files)))
  colnames(circ_counts_df) <- c("sample", "count")
  for (f in 1:length(ciri_files)) {
    sample <- read.table(file = ciri_files[f], 
                         sep = '\t', 
                         skip = 1,
                         header = FALSE)
    
    sample_name <- gsub("\\..*","", ciri_files[f])
    sample_name <- gsub(".*/","", sample_name)
    match_name <- gsub("ccle|gcsi", "",ciri_files[f])
    match_name <- gsub("\\..*","", match_name)
    match_name <- gsub(".*/","", match_name)
    match_name <- gsub("neg","", match_name)
    
    hansen_match <- read.table(file = paste0(RNAseR_dir, match_name,"pos.tsv"), 
                               sep = '\t', 
                               skip = 1,
                               header = FALSE)
    
    colnames(sample) <- c("circRNA_ID","chr","circRNA_start","circRNA_end","junction_reads", "SM_MS_SMS", "non_junction_reads", "junction_reads_ratio", "circRNA_type", "gene_id", "strand","junction_reads_ID")
    colnames(hansen_match) <- c("circRNA_ID","chr","circRNA_start","circRNA_end","junction_reads", "SM_MS_SMS", "non_junction_reads", "junction_reads_ratio", "circRNA_type", "gene_id", "strand","junction_reads_ID")

    sample_match <- sample[which(sample$circRNA_type=="exon"),]
    sample_match <- sample_match[which(sample_match$junction_reads >= 2),]  
    
    hansen_match <- hansen_match[which(hansen_match$circRNA_type=="exon"),]
    hansen_match <- hansen_match[which(hansen_match$junction_reads >= 2),]
    
    sample_match <- sample_match[which(sample_match$circRNA_ID %in% hansen_match$circRNA_ID),]
    hansen_match <- hansen_match[which(hansen_match$circRNA_ID %in% sample_match$circRNA_ID),]
    hansen_match  <- hansen_match[order(sample_match$circRNA_ID),]
    
    filtered <- sample_match[hansen_match$junction_reads/sample_match$junction_reads >= 1.5,]

    count <- sum(filtered$junction_reads)
    circ_counts_df[f,] <- c(sample_name,count)
  }
  return(circ_counts_df)
} 

validated_circRNA_ciri <- filterCIRI(nonRNAseR_dir = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRI2/hansen_match/result", 
                                RNAseR_dir = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRI2/hansen/result/")



validated_circRNA_rminus_ciri <- filterCIRI(nonRNAseR_dir = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRI2/hansen/result", 
                                       RNAseR_dir = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRI2/hansen/result/")


#hansen_matched vs validated hansen_matched bar plot (log2)
library(ggplot2)
cc <- c(rep(c("norm"),4), rep(c("validated"),4))
df.hm <- data.frame("sample"=c(hansen_ciri_matched$sample, validated_circRNA_ciri$sample), "condition"=cc <- c(rep(c("norm"),4), rep(c("validated"),4)), "counts"=log2(as.numeric(c(hansen_ciri_matched$count, validated_circRNA_ciri$count))))

ggplot(df.hm, aes(fill=condition, y=counts, x=sample)) + 
  geom_bar(position="dodge", stat="identity")


#hansen_ribo vs validated hansen_ribo bar plot (log2)
library(ggplot2)
cc <- c(rep(c("norm"),3), rep(c("validated"),3))
df.hm <- data.frame("sample"=c(hansen_ciri_counts[c(1,3,5),"sample"], validated_circRNA_rminus_ciri$sample), "condition"=cc <- c(rep(c("norm"),3), rep(c("validated"),3)), "counts"=log2(as.numeric(c(hansen_ciri_counts[c(1,3,5),"count"], validated_circRNA_ciri$count))))

ggplot(df.hm, aes(fill=condition, y=counts, x=sample)) + 
  geom_bar(position="dodge", stat="identity")





####### read in CIRCexplorer2 abundance of circRNA's #######

#reads in CIRCexplorer2 data
summarizeCIRCexplorer <- function(dir_path){
  
  circ_files <- list.files(path=dir_path,  pattern = ".txt", full.names = TRUE, recursive = TRUE)
  circ_counts_df <- data.frame(matrix(ncol=2, nrow = length(circ_files)))
  colnames(circ_counts_df) <- c("sample", "count")
  for (f in 1:length(circ_files)) {
    sample <- read.table(file = circ_files[f], 
                         sep = '\t',
                         header = FALSE)
    
    colnames(sample) <- c("chrom","start","end", "name","score","strand","thickStart","thickEnd","itemRgb","exonCount","exonSizes", "exonOffsets", "readNumber", "circType", "geneName", "isoformName", "index", "flankIntron")

    sample_name <- sub(".*/ *(.*?) *.txt.*", "\\1", circ_files[f])
    sample <- sample[which(sample$circType=="circRNA"),]
    sample <- sample[which(sample$readNumber >= 2),]
    #count <- length(sample$circType)
    count <- sum(sample$readNumber)
    circ_counts_df[f,] <- c(sample_name,count)
  }
  return(circ_counts_df)
} 

gcsi_circ_counts <- summarizeCIRCexplorer(dir_path = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRCexplorer2/gCSI/unmapped_method/annotate")
ccle_circ_counts <- summarizeCIRCexplorer(dir_path = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRCexplorer2/CCLE/unmapped_method/annotate")
gdsc_circ_counts <- summarizeCIRCexplorer(dir_path = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRCexplorer2/GDSC/unmapped_method/annotate")
hansen_circ_counts <- summarizeCIRCexplorer(dir_path = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRCexplorer2/hansen/annotate")
hansen_circ_matched_counts <- summarizeCIRCexplorer(dir_path = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRCexplorer2/hansen_match/annotate")
#rownames(hansen_ciri_matched)[4] <- "PC3gcsi"

gcsi_circ_counts$sample <- sub(".*gcsi *(.*?) *_.*", "\\1", gcsi_circ_counts$sample)
rownames(gcsi_circ_counts) <- gcsi_circ_counts$sample

ccle_circ_counts$sample <- sub("* *(.*?) *_circular.*", "\\1", ccle_circ_counts$sample)
rownames(ccle_circ_counts) <- ccle_circ_counts$sample

gdsc_circ_counts$sample <- sub(".*gdsc *(.*?) *_circular.*", "\\1", gdsc_circ_counts$sample)
rownames(gdsc_circ_counts) <- gdsc_circ_counts$sample


#average SR technical replicate circRNA count in gCSI
gcsi_circ_counts$cellid <- gcsi$cellid[match(rownames(gcsi_circ_counts), rownames(gcsi))]
gcsi_SR_mean <- mean(as.numeric(gcsi_circ_counts$count[which(gcsi_circ_counts$cellid=="SR")]))
gcsi_circ_counts[which(gcsi_circ_counts$cellid=="SR"),][2,"count"] <- gcsi_SR_mean
gcsi_circ_counts <- gcsi_circ_counts[-which(gcsi_circ_counts$sample=="587641"),]

#order ccle & gdsc cellid by gcsi cellid
ccle_circ_counts$cellid <- ccle$cellid[match(rownames(ccle_circ_counts), rownames(ccle))]
ccle_circ_counts <- ccle_circ_counts[match(gcsi_circ_counts$cellid, ccle_circ_counts$cellid),]

gdsc_circ_counts$cellid <- gdsc$cellid[match(rownames(gdsc_circ_counts), rownames(gdsc))]
gdsc_circ_counts <- gdsc_circ_counts[match(gcsi_circ_counts$cellid, gdsc_circ_counts$cellid),]

circ_combined <- data.frame("gcsi"=as.numeric(gcsi_circ_counts$count), "ccle"=as.numeric(ccle_circ_counts$count), "gdsc"=as.numeric(gdsc_circ_counts$count))
rownames(circ_combined) <- gcsi_circ_counts$cellid




filterCIRC <- function(nonRNAseR_dir, RNAseR_dir){
  #"*neg*" <- for ribominus
  #"\\.txt$"<- for poly-A selected data
  circ_files <- list.files(path=nonRNAseR_dir,  pattern = "*neg*", full.names = T, recursive = TRUE)
  circ_counts_df <- data.frame(matrix(ncol=2, nrow = length(circ_files)))
  colnames(circ_counts_df) <- c("sample", "count")
  for (f in 1:length(circ_files)) {
    sample <- read.table(file = circ_files[f], 
                         sep = '\t',
                         header = FALSE)
    
    sample_name <- gsub("_circularRNA_known.txt","", circ_files[f])
    sample_name <- gsub(".*/","", sample_name)
    match_name <- gsub("neg", "",sample_name)
    #match_name <- gsub("ccle|gcsi", "",sample_name)
    
    hansen_match <- read.table(file = paste0(RNAseR_dir, "/" ,match_name, "pos", "/", match_name,"pos", "_circularRNA_known.txt"),
                               sep = '\t',
                               header = FALSE)
    
    colnames(sample) <- c("chrom","start","end", "name","score","strand","thickStart","thickEnd","itemRgb","exonCount","exonSizes", "exonOffsets", "readNumber", "circType", "geneName", "isoformName", "index", "flankIntron")
    colnames(hansen_match) <- c("chrom","start","end", "name","score","strand","thickStart","thickEnd","itemRgb","exonCount","exonSizes", "exonOffsets", "readNumber", "circType", "geneName", "isoformName", "index", "flankIntron")
    
    sample_match <- sample[which(sample$circType=="circRNA"),]
    sample_match <- sample_match[which(sample_match$readNumber >= 2),]  
    
    hansen_match <- hansen_match[which(hansen_match$circType=="circRNA"),]
    #hansen_match <- hansen_match[which(hansen_match$readNumber >= 2),]
    
    sample_match$pos <- paste0(sample_match$chrom,sample_match$start,sample_match$end)
    hansen_match$pos <- paste0(hansen_match$chrom,hansen_match$start,hansen_match$end)
    
    sample_match <- sample_match[which(sample_match$pos %in% hansen_match$pos),]
    hansen_match <- hansen_match[which(hansen_match$pos %in% sample_match$pos),]
    hansen_match  <- hansen_match[match(sample_match$pos, hansen_match$pos),]
    
    filtered <- sample_match[hansen_match$readNumber/sample_match$readNumber >= 1.5,]
    
    count <- sum(filtered$readNumber)
    circ_counts_df[f,] <- c(sample_name,count)
  }
  return(circ_counts_df)
} 


#tt <- read.table("/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRCexplorer2/hansen_match/annotate/22Rv1ccle/22Rv1ccle_circularRNA_known.txt", sep = '\t', header = FALSE)
#colnames(tt) <- c("chrom","start","end", "name","score","strand","thickStart","thickEnd","itemRgb","exonCount","exonSizes", "exonOffsets", "readNumber", "circType", "geneName", "isoformName", "index", "flankIntron")
#tt$pos <- paste0(tt$chrom,":",tt$start,"-",tt$end)


validated_circRNA_rminus <- filterCIRC(nonRNAseR_dir = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRCexplorer2/hansen/annotate", 
                                       RNAseR_dir = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRCexplorer2/hansen/annotate")

validated_circRNA <- filterCIRC(nonRNAseR_dir = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRCexplorer2/hansen_match/annotate", 
                                       RNAseR_dir = "/Users/anthmam/Desktop/Projects/BHKLAB/ncRNA/results/CIRCexplorer2/hansen/annotate")




####### plot CIRI2 results #######
library(plotly)
data <- data.frame(gcsi_ciri_counts$cellid, ciri_combined$gcsi, ciri_combined$ccle, ciri_combined$gdsc)
data$gcsi_ciri_counts.cellid <- factor(data$gcsi_ciri_counts.cellid, levels = data[["gcsi_ciri_counts.cellid"]])

fig <- plot_ly(data, x = ~gcsi_ciri_counts.cellid, y = ~ciri_combined$gcsi, type = 'bar', name = 'gCSI', marker = list(color = 'rgb(49,130,189)'))
fig <- fig %>% add_trace(y = ~ciri_combined$ccle, name = 'CCLE', marker = list(color = 'rgb(204,204,204)'))
fig <- fig %>% add_trace(y = ~ciri_combined$gdsc, name = 'GDSC', marker = list(color = 'rgb(8,48,107)'))
fig <- fig %>% layout(xaxis = list(title = "", tickangle = -45),
                      yaxis = list(title = ""),
                      margin = list(b = 100),
                      barmode = 'group')


#circRNA RPKM = junction reads/(circRNA length × total mapped reads)

