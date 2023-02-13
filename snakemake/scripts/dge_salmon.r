### (c) Shweta V. Pipaliya 
### 2022-04-09 
## This script uses DESEQ2 to perform DGE using gene_level counts generated either using Salmon
### How to run from command-line: Rscript --vanilla DGE_Salmon_DESEQ2.svp.R
### Bioconductor packages in this script require R v.4.2.0

## Set seed 
set.seed(777)

## Set current working directory to output DGE summary results and plots
wd <- getwd() 
if (!is.null(wd)) setwd(wd) # set current directory as working directory; wd should have the two input files and this rscript
dir.create("dge_salmon.results.svp") # create a new results folder 
list.files(wd)

## Install and load required libraries (addition of tximport, readr, and biomaRt)
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

#BiocManager::install("readr")
library(readr)

#BiocManager::install("biomaRt")
library(biomaRt)

#BiocManager::install("tibble")
library(tibble)

#BiocManager::install("apeglm")
library(apeglm)

#BiocManager::install("DESeq2")
library(DESeq2)

#BiocManager::install("BiocParallel")
library(BiocParallel)
register(MulticoreParam(4)) # register the number of nodes on multicoreParam for DESeq2 run

#if (!require(pheatmap)) install.packages('pheatmap')
library(pheatmap)

#if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)

#if (!require(ggplot2)) install.packages('ggplot2')
library(ggplot2)

#if (!require(ggrepel)) install.packages('ggrepel')
library(ggrepel)

#if (!require(RColorBrewer)) install.packages('RColorBrewer')
library(RColorBrewer)

#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

#BiocManager::install("DOSE")
library(DOSE)

#BiocManager::install("enrichplot")
library(enrichplot)

organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE, force = TRUE)
library(organism, character.only = TRUE)

#install.packages("ggnewscale")
library(ggnewscale)

#install.packages("ggupset")
library(ggupset)

#BiocManager::install("pathview")
library(pathview)

#BiocManager::install("clusterProfiler")
library(clusterProfiler)

#BiocManager::install("AnnotationHub")
library(AnnotationHub)

#BiocManager::install("ensembldb")
library(ensembldb)

install.packages("devtools")
devtools::install_github("stephenturner/annotables")

### FIRST PART OF THE SCRIPT: Perform sample and gene qc, count normalization, and DGE using DESeq2 ###
## load phenotypes with samples corresponding to respective disease groups (asymptomatic, mild, severe, critical)
metaData <- read.csv("../PRJEB43380.metadata.csv", header=TRUE, sep=",", row.names = 1)
#mt_df1 <- data.frame(metaData[,-1], row.names= metaData[,1]) # in case there are duplicate rows then use this command to determine which one

## check the sample names
colnames(txi$counts)

## create object of class DESeqDataSet to create a matrix that can be read by DESeq2
deseq2Data <- DESeqDataSetFromTximport(txi, 
                                colData = metaData, 
                                design = ~ sample_group)
print(deseq2Data)

## transform gene counts for PCA for sample-level QC 
pdf("dge_salmon.results.svp/plot.PRJEB43380.sample.pca.pdf")

rld <- rlog(deseq2Data, blind=TRUE)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat)) # run prcomp for PCA
df.pca = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], 
                    PC3 = pca$x[,3], PC4 = pca$x[,4]) 

plot.pca <- plotPCA(rld,intgroup="sample_group") # use DESeq2's plotPCA() for plotting
plot.pca + ggtitle("Principal Component Analysis: Visualizing sample variance")

dev.off()

## Remove sample outliers - perform this step only if necessary
outliers <- as.character(subset(colnames(deseq2Data), df.pca$PC1 > 100))
print (outliers)

deseq2Data.sqc <- deseq2Data[, !(colnames(deseq2Data) %in% outliers)] 

## extract the rlog matrix from the QC'd deseq2 data object and re-run pca + plot to see how the results looks now
pdf("dge_salmon.results.svp/plot.PRJEB43380.qc.sample.pca.pdf")

rld.sqc <- rlog(deseq2Data.sqc, blind=T)
rld_mat.sqc <- assay(rld.sqc)   

plot.pca <- plotPCA(rld.sqc,intgroup="sample_group") # use DESeq2's plotPCA() for plotting
plot.pca + ggtitle("Principal Component Analysis: Visualizing sample variance")

dev.off()

## compute pairwise correlation values for post-QC samples and plot a heatmap
pdf("dge_salmon.results.svp/plot.sampleQC.heatmap.pdf")

rld_cor <- cor(rld_mat.sqc)   
head(rld_cor)

heat.colors <- brewer.pal(6, "Blues")
pheatmap(rld_cor, 
         main = "Heatmap: Pairwise correlation between samples",
         color = heat.colors, 
         border_color=NA, 
         fontsize = 5, 
         fontsize_row = 5,
         height=20)

dev.off()

## perform gene-level QC by filtering out rows with less than 10 reads across all the samples
dim(deseq2Data)
dim(deseq2Data[rowSums(counts(deseq2Data)) > 10, ]) # assess filter effect
deseq2Data.gqc <- deseq2Data[rowSums(counts(deseq2Data)) > 10, ] # run and output filtered data

## perform DGE analyses with pairwise comparisons on the filtered dataset using the DESeq() on multiple cores - may take several minutes

deseq2Data_DGE <- DESeq(deseq2Data.gqc, parallel = TRUE)

# DGE for Asymptomatic vs mild
AM <- results(deseq2Data_DGE, contrast=c("sample_group", "asymptomatic", "mild"), 
              independentFiltering=TRUE, alpha = 0.05, 
              pAdjustMethod="BH")
#AM <- lfcShrink(deseq2Data_DGE, contrast=c("sample_group", "asymptomatic", "mild"), coef="sample_group_asymptomatic_vs_mild", type="apeglm", res=AM)
dgeSum_AM <- summary(AM) # view summary summary statistics for up- and down-regulated genes
head(AM)

dgeResults_AM <- AM[order(AM$padj),] ## sort and write out DGE results for all genes by adjusted p-values (Benjamini-Hochberg FDR method)
head(dgeResults_AM)

write.csv(as.data.frame(AM[order(AM$padj),]), file = "dge_salmon.results.svp/asymptomatic_v_mild_dge.csv")

# DGE for asymptomatic vs severe
AS <- results(deseq2Data_DGE, contrast=c("sample_group", "asymptomatic", "severe"), 
              independentFiltering=TRUE, alpha = 0.05, 
              pAdjustMethod="BH")
#AS <- lfcShrink(deseq2Data_DGE, contrast=c("sample_group", "asymptomatic", "severe"), res=AS)
dgeSum_AS <- summary(AS) # view summary summary statistics for up- and down-regulated genes
head(AS)

dgeResults_AS <- AS[order(AS$padj),] ## sort and write out DGE results for all genes by adjusted p-values (Benjamini-Hochberg FDR method)
head(dgeResults_AS)

write.csv(as.data.frame(AS[order(AS$padj),]), file = "dge_salmon.results.svp/asymptomatic_v_severe_dge.csv")

# DGE for asymptomatic vs critical
AC <- results(deseq2Data_DGE, contrast=c("sample_group", "asymptomatic", "critical"), 
              independentFiltering=TRUE, alpha = 0.05, 
              pAdjustMethod="BH")
#AC <- lfcShrink(deseq2Data_DGE, contrast=c("sample_group", "asymptomatic", "critical"), res=AC)
dgeSum_AC <- summary(AC) # view summary summary statistics for up- and down-regulated genes
head(AC)

dgeResults_AC <- AC[order(AC$padj),] ## sort and write out DGE results for all genes by adjusted p-values (Benjamini-Hochberg FDR method)
head(dgeResults_AC)

write.csv(as.data.frame(AC[order(AC$padj),]), file = "dge_salmon.results.svp/asymptomatic_v_critical_dge.csv")

# DGE for mild vs severe
MS <- results(deseq2Data_DGE, contrast=c("sample_group", "mild", "severe"), 
              independentFiltering=TRUE, alpha = 0.05, 
              pAdjustMethod="BH")
#MS <- lfcShrink(deseq2Data_DGE, contrast=c("sample_group", "mild", "severel"), res=MS)
dgeSum_MS <- summary(MS) # view summary summary statistics for up- and down-regulated genes
head(MS)

dgeResults_MS <- MS[order(MS$padj),] ## sort and write out DGE results for all genes by adjusted p-values (Benjamini-Hochberg FDR method)
head(dgeResults_MS)

write.csv(as.data.frame(MS[order(MS$padj),]), file = "dge_salmon.results.svp/mild_v_severe_dge.csv")

# DGE for mild vs critical
MC <- results(deseq2Data_DGE, contrast=c("sample_group", "mild", "critical"), 
              independentFiltering=TRUE, alpha = 0.05, 
              pAdjustMethod="BH")
#MC <- lfcShrink(deseq2Data_DGE, contrast=c("sample_group", "mild", "critical"), res=MC)
dgeSum_MC <- summary(MC) # view summary summary statistics for up- and down-regulated genes
head(MC)

dgeResults_MC <- MC[order(MC$padj),] ## sort and write out DGE results for all genes by adjusted p-values (Benjamini-Hochberg FDR method)
head(dgeResults_MC)

write.csv(as.data.frame(MC[order(MC$padj),]), file = "dge_salmon.results.svp/mild_v_critical_dge.csv")

# DGE for severe vs critical
SC <- results(deseq2Data_DGE, contrast=c("sample_group", "severe", "critical"), 
              independentFiltering=TRUE, alpha = 0.05, 
              pAdjustMethod="BH")
#SC <- lfcShrink(deseq2Data_DGE, contrast=c("sample_group", "severe", "critical"), res=SC)
dgeSum_SC <- summary(SC) # view summary summary statistics for up- and down-regulated genes
head(SC)

dgeResults_SC <- SC[order(SC$padj),] ## sort and write out DGE results for all genes by adjusted p-values (Benjamini-Hochberg FDR method)
head(dgeResults_SC)

write.csv(as.data.frame(SC[order(SC$padj),]), file = "dge_salmon.results.svp/severe_v_critical_dge.csv")

## get gene names for the top 50 dge for AS
top50_gene_AS <- rownames(dgeResults_AS[1:50, ]) # get rownames
top50_sigOE_AS_dge <- dgeResults_AS[1:50, ] # extract dge results for the top 50 genes for AS

write.csv(as.data.frame(dgeResults_AS[1:50, ] ), file="dge_salmon.results.svp/top50.AS_dge.csv")

## get gene names for the top 50 dge for AM
top50_gene_AM <- rownames(dgeResults_AM[1:50, ]) # get rownames
top50_sigOE_AM_dge <- dgeResults_AM[1:50, ] # extract dge results for the top 50 genes for AM

write.csv(as.data.frame(dgeResults_AM[1:50, ] ), file="dge_salmon.results.svp/top50.AM_dge.csv")

## get gene names for the top 50 dge for AC
top50_gene_AC <- rownames(dgeResults_AC[1:50, ]) # get rownames
top50_sigOE_AC_dge <- dgeResults_AC[1:50, ] # extract dge results for the top 50 genes for AC

write.csv(as.data.frame(dgeResults_AC[1:50, ] ), file="dge_salmon.results.svp/top50.AC_dge.csv")

## get gene names for the top 50 dge for MS
top50_gene_MS <- rownames(dgeResults_MS[1:50, ]) # get rownames
top50_sigOE_MS_dge <- dgeResults_MS[1:50, ] # extract dge results for the top 50 genes for MS

write.csv(as.data.frame(dgeResults_MS[1:50, ] ), file="dge_salmon.results.svp/top50.MS_dge.csv")

## get gene names for the top 50 dge for MC
top50_gene_MC <- rownames(dgeResults_MC[1:50, ]) # get rownames
top50_sigOE_MC_dge <- dgeResults_MC[1:50, ] # extract dge results for the top 50 genes for MC

write.csv(as.data.frame(dgeResults_MC[1:50, ] ), file="dge_salmon.results.svp/top50.MC_dge.csv")

## get gene names for the top 50 dge for SC
top50_gene_SC <- rownames(dgeResults_SC[1:50, ]) # get rownames
top50_sigOE_SC_dge <- dgeResults_SC[1:50, ] # extract dge results for the top 50 genes for SC

write.csv(as.data.frame(dgeResults_SC[1:50, ] ), file="dge_salmon.results.svp/top50.SC_dge.csv")

## get and write out normalized counts for all genes in all conditions
normalized_counts <- counts(deseq2Data_DGE, normalized=TRUE)
head(normalized_counts)

write.csv(as.data.frame(counts(deseq2Data_DGE, normalized=TRUE) ), file="dge_salmon.results.svp/all.normalized_counts.csv")

## get and write out normalized counts for the top 50 differentially expressed genes for AS
top50_sigOE_AS_norm <- normalized_counts[top50_gene_AS, ]
write.csv(as.data.frame(normalized_counts[top50_gene_AS, ] ), file="dge_salmon.results.svp/top50.AS_nc.csv")

## get and write out normalized counts for the top 50 differentially expressed genes for AM
top50_sigOE_AM_norm <- normalized_counts[top50_gene_AM, ]
write.csv(as.data.frame(normalized_counts[top50_gene_AM, ] ), file="dge_salmon.results.svp/top50.AM_nc.csv")

## get and write out normalized counts for the top 50 differentially expressed genes for AC
top50_sigOE_AC_norm <- normalized_counts[top50_gene_AC, ]
write.csv(as.data.frame(normalized_counts[top50_gene_AC, ] ), file="dge_salmon.results.svp/top50.AC_nc.csv")

## get and write out normalized counts for the top 50 differentially expressed genes for MS
top50_sigOE_MS_norm <- normalized_counts[top50_gene_MS, ]
write.csv(as.data.frame(normalized_counts[top50_gene_MS, ] ), file="dge_salmon.results.svp/top50.MS_nc.csv")

## get and write out normalized counts for the top 50 differentially expressed genes for MC
top50_sigOE_MC_norm <- normalized_counts[top50_gene_MC, ]
write.csv(as.data.frame(normalized_counts[top50_gene_MC, ] ), file="dge_salmon.results.svp/top50.MC_nc.csv")

## get and write out normalized counts for the top 50 differentially expressed genes for SC
top50_sigOE_SC_norm <- normalized_counts[top50_gene_SC, ]
write.csv(as.data.frame(normalized_counts[top50_gene_SC, ] ), file="dge_salmon.results.svp/top50.SC_nc.csv")

### SECOND PART OF THE SCRIPT: Visualize log fold change and normalized counts for biomarkers of interest ###

## generate an MA plot for all p-adj DGE results for AM
pdf(file = "dge_salmon.results.svp/plot.MA.AM.pdf") 

df.AM.deSeqRes <- data.frame(dgeResults_AM, stringsAsFactors = FALSE)
df.AM.deSeqRes <- tibble::rownames_to_column(df.AM.deSeqRes, "row_names") # convert original dgeResults df with gene row names as a separate column for MA plotting
colnames(df.AM.deSeqRes)[1] <- "Gene" # change column 1name to Gene

df.AM.deSeqRes$Significance <- ifelse(df.AM.deSeqRes$padj <= 0.05, "Yes", "No") # specify if gene significance meets padj =< 0.05 threshold
df.AM.deSeqRes$genelabels <- factor(df.AM.deSeqRes$Gene, levels = c(top50_gene_AM)) # add label column 

ggplot(df.AM.deSeqRes) +
  geom_point(aes(x = log10(baseMean), 
                 y = log2FoldChange, 
                 color = Significance), 
             size = 0.15) +
  geom_label_repel(aes(x = log10(baseMean), 
                       y = log2FoldChange,
                       label = genelabels),
                   nudge_x = 0.1,
                   nudge_y = 0.45, 
                   segment.curvature = -1e-20,
                   size = 1) +
  xlab("log10(baseMean)") +
  ylab("log2FoldChange") +
  ylim(-10,10) +
  geom_hline(yintercept = c(0.58, -0.58),       #optionally add a fold-change threshold line
             linetype = "dashed") +
  ggtitle("MA plot AM: Log fold-change vs. mean expression") +
  theme(plot.title = element_text(size = rel(1.0), hjust = 0.5),
        axis.title = element_text(size = rel(0.75))) 

dev.off()

## generate an MA plot for all p-adj DGE results for AC
pdf(file = "dge_salmon.results.svp/plot.MA.AC.pdf") 

df.AC.deSeqRes <- data.frame(dgeResults_AC, stringsAsFactors = FALSE)
df.AC.deSeqRes <- tibble::rownames_to_column(df.AC.deSeqRes, "row_names") # convert original dgeResults df with gene row names as a separate column for MA plotting
colnames(df.AC.deSeqRes)[1] <- "Gene" # change column 1name to Gene

df.AC.deSeqRes$Significance <- ifelse(df.AC.deSeqRes$padj <= 0.05, "Yes", "No") # specify if gene significance meets padj =< 0.05 threshold
df.AC.deSeqRes$genelabels <- factor(df.AC.deSeqRes$Gene, levels = c(top50_gene_AC)) # add label column 

ggplot(df.AC.deSeqRes) +
  geom_point(aes(x = log10(baseMean), 
                 y = log2FoldChange, 
                 color = Significance), 
             size = 0.15) +
  geom_label_repel(aes(x = log10(baseMean), 
                       y = log2FoldChange,
                       label = genelabels),
                   nudge_x = 0.1,
                   nudge_y = 0.45, 
                   segment.curvature = -1e-20,
                   size = 1) +
  xlab("log10(baseMean)") +
  ylab("log2FoldChange") +
  ylim(-10,10) +
  geom_hline(yintercept = c(0.58, -0.58),       #optionally add a fold-change threshold line
             linetype = "dashed") +
  ggtitle("MA plot AC: Log fold-change vs. mean expression") +
  theme(plot.title = element_text(size = rel(1.0), hjust = 0.5),
        axis.title = element_text(size = rel(0.75))) 

dev.off()

## generate an MA plot for all p-adj DGE results for AS
pdf(file = "dge_salmon.results.svp/plot.MA.AS.pdf") 

df.AS.deSeqRes <- data.frame(dgeResults_AS, stringsAsFactors = FALSE)
df.AS.deSeqRes <- tibble::rownames_to_column(df.AS.deSeqRes, "row_names") # convert original dgeResults df with gene row names as a separate column for MA plotting
colnames(df.AS.deSeqRes)[1] <- "Gene" # change column 1name to Gene

df.AS.deSeqRes$Significance <- ifelse(df.AS.deSeqRes$padj <= 0.05, "Yes", "No") # specify if gene significance meets padj =< 0.05 threshold
df.AS.deSeqRes$genelabels <- factor(df.AS.deSeqRes$Gene, levels = c(top50_gene_AS)) # add label column 

ggplot(df.AS.deSeqRes) +
  geom_point(aes(x = log10(baseMean), 
                 y = log2FoldChange, 
                 color = Significance), 
             size = 0.15) +
  geom_label_repel(aes(x = log10(baseMean), 
                       y = log2FoldChange,
                       label = genelabels),
                   nudge_x = 0.1,
                   nudge_y = 0.45, 
                   segment.curvature = -1e-20,
                   size = 1) +
  xlab("log10(baseMean)") +
  ylab("log2FoldChange") +
  ylim(-10,10) +
  geom_hline(yintercept = c(0.58, -0.58),       #optionally add a fold-change threshold line
             linetype = "dashed") +
  ggtitle("MA plot AS: Log fold-change vs. mean expression") +
  theme(plot.title = element_text(size = rel(1.0), hjust = 0.5),
        axis.title = element_text(size = rel(0.75))) 

dev.off()

## generate an MA plot for all p-adj DGE results for MS
pdf(file = "dge_salmon.results.svp/plot.MA.MS.pdf") 

df.MS.deSeqRes <- data.frame(dgeResults_MS, stringsAsFactors = FALSE)
df.MS.deSeqRes <- tibble::rownames_to_column(df.MS.deSeqRes, "row_names") # convert original dgeResults df with gene row names as a separate column for MA plotting
colnames(df.MS.deSeqRes)[1] <- "Gene" # change column 1name to Gene

df.MS.deSeqRes$Significance <- ifelse(df.MS.deSeqRes$padj <= 0.05, "Yes", "No") # specify if gene significance meets padj =< 0.05 threshold
df.MS.deSeqRes$genelabels <- factor(df.MS.deSeqRes$Gene, levels = c(top50_gene_MS)) # add label column 

ggplot(df.MS.deSeqRes) +
  geom_point(aes(x = log10(baseMean), 
                 y = log2FoldChange, 
                 color = Significance), 
             size = 0.15) +
  geom_label_repel(aes(x = log10(baseMean), 
                       y = log2FoldChange,
                       label = genelabels),
                   nudge_x = 0.1,
                   nudge_y = 0.45, 
                   segment.curvature = -1e-20,
                   size = 1) +
  xlab("log10(baseMean)") +
  ylab("log2FoldChange") +
  ylim(-10,10) +
  geom_hline(yintercept = c(0.58, -0.58),       #optionally add a fold-change threshold line
             linetype = "dashed") +
  ggtitle("MA plot MS: Log fold-change vs. mean expression") +
  theme(plot.title = element_text(size = rel(1.0), hjust = 0.5),
        axis.title = element_text(size = rel(0.75))) 

dev.off()

## generate an MA plot for all p-adj DGE results for MC
pdf(file = "dge_salmon.results.svp/plot.MA.MC.pdf") 

df.MC.deSeqRes <- data.frame(dgeResults_MC, stringsAsFactors = FALSE)
df.MC.deSeqRes <- tibble::rownames_to_column(df.MC.deSeqRes, "row_names") # convert original dgeResults df with gene row names as a separate column for MA plotting
colnames(df.MC.deSeqRes)[1] <- "Gene" # change column 1name to Gene

df.MC.deSeqRes$Significance <- ifelse(df.MC.deSeqRes$padj <= 0.05, "Yes", "No") # specify if gene significance meets padj =< 0.05 threshold
df.MC.deSeqRes$genelabels <- factor(df.MC.deSeqRes$Gene, levels = c(top50_gene_MC)) # add label column 

ggplot(df.MC.deSeqRes) +
  geom_point(aes(x = log10(baseMean), 
                 y = log2FoldChange, 
                 color = Significance), 
             size = 0.15) +
  geom_label_repel(aes(x = log10(baseMean), 
                       y = log2FoldChange,
                       label = genelabels),
                   nudge_x = 0.1,
                   nudge_y = 0.45, 
                   segment.curvature = -1e-20,
                   size = 1) +
  xlab("log10(baseMean)") +
  ylab("log2FoldChange") +
  ylim(-10,10) +
  geom_hline(yintercept = c(0.58, -0.58),       #optionally add a fold-change threshold line
             linetype = "dashed") +
  ggtitle("MA plot MC: Log fold-change vs. mean expression") +
  theme(plot.title = element_text(size = rel(1.0), hjust = 0.5),
        axis.title = element_text(size = rel(0.75))) 

dev.off()

## generate an MA plot for all p-adj DGE results for SC
pdf(file = "dge_salmon.results.svp/plot.MA.SC.pdf") 

df.SC.deSeqRes <- data.frame(dgeResults_SC, stringsAsFactors = FALSE)
df.SC.deSeqRes <- tibble::rownames_to_column(df.SC.deSeqRes, "row_names") # convert original dgeResults df with gene row names as a separate column for MA plotting
colnames(df.SC.deSeqRes)[1] <- "Gene" # change column 1name to Gene

df.SC.deSeqRes$Significance <- ifelse(df.SC.deSeqRes$padj <= 0.05, "Yes", "No") # specify if gene significance meets padj =< 0.05 threshold
df.SC.deSeqRes$genelabels <- factor(df.SC.deSeqRes$Gene, levels = c(top50_gene_SC)) # add label column 

ggplot(df.SC.deSeqRes) +
  geom_point(aes(x = log10(baseMean), 
                 y = log2FoldChange, 
                 color = Significance), 
             size = 0.15) +
  geom_label_repel(aes(x = log10(baseMean), 
                       y = log2FoldChange,
                       label = genelabels),
                   nudge_x = 0.1,
                   nudge_y = 0.45, 
                   segment.curvature = -1e-20,
                   size = 1) +
  xlab("log10(baseMean)") +
  ylab("log2FoldChange") +
  ylim(-10,10) +
  geom_hline(yintercept = c(0.58, -0.58),       #optionally add a fold-change threshold line
             linetype = "dashed") +
  ggtitle("MA plot SC: Log fold-change vs. mean expression") +
  theme(plot.title = element_text(size = rel(1.0), hjust = 0.5),
        axis.title = element_text(size = rel(0.75))) 

dev.off()

## Generate a volcano plot using the EnhancedVolcano plot package for AM
dge.AM.vp <- read.csv("dge_salmon.results.svp/asymptomatic_v_mild_dge.csv", row.names = 1, sep=",")

pdf(file = "dge_salmon.results.svp/plot.volcano.AM.pdf") 

EnhancedVolcano(dge.AM.vp,
                lab = rownames(dge.AM.vp),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "Volcano Plot",
                subtitle = "Differential expression for asymptomatic v. mild COVID-19", 
                pCutoff = 1.34e-06,
                FCcutoff = 1.5)
dev.off()

## Generate a volcano plot using the EnhancedVolcano plot package for AS
dge.AS.vp <- read.csv("dge_salmon.results.svp/asymptomatic_v_severe_dge.csv", row.names = 1, sep=",")

pdf(file = "dge_salmon.results.svp/plot.volcano.AS.pdf") 

EnhancedVolcano(dge.AS.vp,
                lab = rownames(dge.AS.vp),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "Volcano Plot",
                subtitle = "Differential expression for asymptomatic v. severe COVID-19", 
                pCutoff = 1.34e-06,
                FCcutoff = 1.5)

dev.off()

## Generate a volcano plot using the EnhancedVolcano plot package for AC
dge.AC.vp <- read.csv("dge_salmon.results.svp/asymptomatic_v_critical_dge.csv", row.names = 1, sep=",")

pdf(file = "dge_salmon.results.svp/plot.volcano.AC.pdf") 

EnhancedVolcano(dge.AC.vp,
                lab = rownames(dge.AC.vp),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "Volcano Plot",
                subtitle = "Differential expression for asymptomatic v. critical COVID-19", 
                pCutoff = 1.34e-06,
                FCcutoff = 1.5)

dev.off()

## Generate a volcano plot using the EnhancedVolcano plot package for MS
dge.MS.vp <- read.csv("dge_salmon.results.svp/mild_v_severe_dge.csv", row.names = 1, sep=",")

pdf(file = "dge_salmon.results.svp/plot.volcano.MS.pdf") 

EnhancedVolcano(dge.MS.vp,
                lab = rownames(dge.MS.vp),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "Volcano Plot",
                subtitle = "Differential expression for mild v. severe COVID-19", 
                pCutoff = 1.34e-06,
                FCcutoff = 1.5)

dev.off()

## Generate a volcano plot using the EnhancedVolcano plot package for MC
dge.MC.vp <- read.csv("dge_salmon.results.svp/mild_v_critical_dge.csv", row.names = 1, sep=",")

pdf(file = "dge_salmon.results.svp/plot.volcano.MC.pdf") 

EnhancedVolcano(dge.MC.vp,
                lab = rownames(dge.MC.vp),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "Volcano Plot",
                subtitle = "Differential expression for mild v. critical COVID-19", 
                pCutoff = 1.34e-06,
                FCcutoff = 1.5)

dev.off()

## Generate a volcano plot using the EnhancedVolcano plot package for SC
dge.SC.vp <- read.csv("dge_salmon.results.svp/severe_v_critical_dge.csv", row.names = 1, sep=",")

pdf(file = "dge_salmon.results.svp/plot.volcano.SC.pdf") 

EnhancedVolcano(dge.SC.vp,
                lab = rownames(dge.SC.vp),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "Volcano Plot",
                subtitle = "Differential expression for severe v. critical COVID-19", 
                pCutoff = 1.34e-06,
                FCcutoff = 1.5)

dev.off()

## plot boxplots of normalized count comparisons for top 20 biomarkers in healthy vs. disease samples
pdf(file = "dge_salmon.results.svp/plot.boxplot.pdf")

goi_NC <- read.csv("/Users/ShwetaPipaliya/Documents/PostDoc/projects/viral_transcriptomics/analyses/sars_cov_2/GSE172114/06_DGE/dge.results.svp/top20.nc_v_critical_nc.csv") # read in the normalized counts csv for biomarkers from previous steps
metaData_plot <- read.csv("/Users/ShwetaPipaliya/Documents/PostDoc/projects/viral_transcriptomics/analyses/sars_cov_2/GSE172114/06_DGE/metadata.deseq.csv") # read in the metadata for table joining 

colnames(goi_NC)[1] <- "gene" # change column 1 name to gene for plotting
gathered_goi_NC <- goi_NC %>% # gather normalized counts for each sample (cols 2-6) into a single column
  gather(colnames(goi_NC) [2:50], key = "sample", value = "normalized_counts")

gathered_goi_NC <- inner_join(metaData_plot, gathered_goi_NC)  # join metadata (healthy vs. disease) with the gathered df

ggplot(gathered_goi_NC, 
       aes(x = gene, 
           y = normalized_counts)) +
  scale_y_log10() +
  ylab("Normalized Counts") +
  xlab("Gene") +
  ggtitle("Box plot: Normalized counts for top 20 DGE genes") + 
  geom_boxplot(aes(fill = sample_group), 
               position = position_dodge(0.9)) +
  facet_wrap(~ gene, scales = "free")

dev.off()

## draw a heatmap for Z-scores of the log-scaled normalized gene counts for the biomarkers - Not Plotted Yet
pdf(file="/Users/ShwetaPipaliya/Documents/PostDoc/projects/viral_transcriptomics/analyses/sars_cov_2/GSE172114/06_DGE/dge.results.svp/plot.nc.heatmap.pdf")

heatmap_data <- goi_NC %>%
  select(1:50) %>%
  column_to_rownames("gene") 

annotation <- metaData_plot %>%    # set annotation for the individual sample names "Disease" and "Healthy"
  select(sample, sample_group) %>%
  data.frame(row.names = "sample")

pheatmap(log2 (heatmap_data + 1), 
         main = "Heatmap: Z-scores of log transformed normalized counts",
         cluster_rows = T, 
         color = hcl.colors(50, "BluYl"),
         show_rownames = T, 
         annotation = annotation,
         display_numbers = FALSE,
         number_color = "black",
         fontsize_number = 5, 
         border_color = "grey", 
         fontsize = 5, 
         cluster_cols = FALSE,
         scale = "row",          # "scale = row" parameter computes z-scores AFTER clustering
         fontsize_row = 5,
         height = 10)

dev.off()

### END ###
