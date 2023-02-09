### (c) Shweta V. Pipaliya 
### 2022-04-09 
## This script performs DGE analyses using Salmon quant.sf data using DESEQ2
## Overall workflow based on HBCtraining tutorial: https://hbctraining.github.io/DGE_workshop_salmon/lessons/01_DGE_setup_and_overview.html
### How to run this script from command-line: Rscript --vanilla 10_DGE_Salmon_DESEQ2.svp.R
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

#BiocManager::install("tximport")
library(tximport)

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

### FIRST PART OF THE SCRIPT: Format the Salmon TPM to gene-level counts using tximport, sample and gene qc, count normalization, and DGE using DESeq2 ###
## Convert salmon TPM to gene-level counts using tximport

## List and name all the quant.sf file for each sample  
samples <- Sys.glob("*.quant.sf")
names(samples) <- samples
print(samples)
all(file.exists(samples))

## Create a character vector of Ensembl IDs
ids <- read.delim(samples[1], sep="\t", header=T) # extract the transcript ids from one of the files
ids <- as.character(ids[,1])
require(stringr)
ids.strip <- str_replace(ids, "([.][0-9])", "") # this command removes version numbers from each accession

## Obtain a vector of the combined salmon quant text file to extract the gene names - NULL
#salmon.quant <- read.delim("/Users/ShwetaPipaliya/Documents/PostDoc/projects/viral_transcriptomics/analyses/sars_cov_2/PRJEB43380/salmon_quant/salmon_combined_counts.txt", header=TRUE, sep="\t")
#sq <- tibble::rownames_to_column(salmon.quant, "refseq_id")

## Convert Refseq IDs to geneIds using biomaRt and the ENSEMBL dataset (this step may take several minutes)
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL", host="https://mar2016.archive.ensembl.org")) # in case the below command times out at the ENSEMBL mirror
mart <- biomaRt::useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')
tx2gene <- getBM(
  filters="refseq_mrna", 
  attributes = c("refseq_mrna","hgnc_symbol"),
  values=ids.strip, 
  mart= mart)

tx2gene %>% View() # take a look at the resulting tx2gene output

## Run tximport to summarize gene-level information
?tximport   # take a look at the arguments for the tximport function
txi <- tximport(samples, type="salmon", tx2gene=tx2gene[,c("refseq_mrna", "hgnc_symbol")], countsFromAbundance="lengthScaledTPM", ignoreTxVersion=TRUE)
attributes(txi) # read the output from the txi object 

head(txi$counts)

## Perform differential gene expression analyses using DESEQ2
## load sample metadata with samples corresponding to respective disease groups (asymptomatic, mild, severe, critical)
metaData <- read.csv("/Users/ShwetaPipaliya/Documents/PostDoc/projects/viral_transcriptomics/analyses/sars_cov_2/PRJEB43380/DGE_salmon/PRJEB43380.metadata.csv", header=TRUE, sep=",", row.names = 1)
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

## Remove sample outliers clustering away - perform this step only if necessary
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

## plot boxplots of normalized count comparisons for top 20 biomarkers in healthy vs. disease samples - NOT PLOTTED YET
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

## Perform OverExpression Analysis (ORA) using Clusterprofiler and Gene Ontology

## clusterProfiler does not work as easily using gene names, so we will turn gene names into Ensembl IDs using 
## clusterProfiler::bitr and merge the IDs back with the DE results
keytypes(org.Hs.eg.db)

## do this for the AM comparision
ids_AM <- bitr(rownames(dge.AM.vp), 
            fromType = "SYMBOL", 
            toType = c("ENSEMBL", "ENTREZID"), 
            OrgDb = "org.Hs.eg.db")

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates_AM <- which(duplicated(ids_AM$SYMBOL) == FALSE)

ids_AM <- ids_AM[non_duplicates_AM, ] 

## Merge the Ensembl IDs with the results     
merged_gene_ids_AM <- merge(x=dge.AM.vp, y=ids_AM, by.x="row.names", by.y="SYMBOL")             

## Extract significant results
AM_sigOE <- subset(merged_gene_ids_AM, padj < 0.05)

AM_sigOE_genes <- as.character(AM_sigOE$ENSEMBL)

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
AM_allOE_genes <- as.character(merged_gene_ids_AM$ENSEMBL)

## Run GO enrichment analysis 
ego_AM <- enrichGO(gene = AM_sigOE_genes, 
                universe = AM_allOE_genes, 
                keyType = "ENSEMBL", 
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

## Output results from GO analysis to a table
cluster_summary_AM <- data.frame(ego_AM)

write.csv(cluster_summary_AM, "/Users/ShwetaPipaliya/Documents/PostDoc/projects/viral_transcriptomics/analyses/sars_cov_2/PRJEB43380/enrichment_analyses/clusterProfiler_AM.csv")

## Dotplot gives top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not padj.
dotplot(ego_AM, showCategory=50)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
AM_OE_foldchanges <- AM_sigOE$log2FoldChange

names(AM_OE_foldchanges) <- AM_sigOE$Row.names

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego_AM, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=AM_OE_foldchanges, 
         vertex.label.font=6)

## If some of the high fold changes are getting drowned out due to a large range, you could set a minimum and maximum fold change value
AM_OE_foldchanges <- ifelse(AM_OE_foldchanges > 2, 2, AM_OE_foldchanges)
AM_OE_foldchanges <- ifelse(AM_OE_foldchanges < -2, -2, AM_OE_foldchanges)

cnetplot(ego_AM, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=AM_OE_foldchanges, 
         vertex.label.font=6)

## Repeat for the AS comparision
ids_AS <- bitr(rownames(dge.AS.vp), 
            fromType = "SYMBOL", 
            toType = c("ENSEMBL", "ENTREZID"), 
            OrgDb = "org.Hs.eg.db")

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates_AS <- which(duplicated(ids_AS$SYMBOL) == FALSE)

ids_AS <- ids_AS[non_duplicates_AS, ] 

## Merge the Ensembl IDs with the results     
merged_gene_ids_AS <- merge(x=dge.AS.vp, y=ids_AS, by.x="row.names", by.y="SYMBOL")             

## Extract significant results
AS_sigOE <- subset(merged_gene_ids_AS, padj < 0.05)

AS_sigOE_genes <- as.character(AS_sigOE$ENSEMBL)

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
AS_allOE_genes <- as.character(merged_gene_ids_AS$ENSEMBL)

## Run GO enrichment analysis 
ego_AS <- enrichGO(gene = AS_sigOE_genes, 
                universe = AS_allOE_genes, 
                keyType = "ENSEMBL", 
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

## Output results from GO analysis to a table
cluster_summary_AS <- data.frame(ego_AS)

write.csv(cluster_summary_AS, "/Users/ShwetaPipaliya/Documents/PostDoc/projects/viral_transcriptomics/analyses/sars_cov_2/PRJEB43380/enrichment_analyses/clusterProfiler_AS.csv")

## Dotplot gives top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not padj.
dotplot(ego_AS, showCategory=50)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
AS_OE_foldchanges <- AS_sigOE$log2FoldChange

names(AS_OE_foldchanges) <- AS_sigOE$Row.names

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego_AS, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=AS_OE_foldchanges, 
         vertex.label.font=6)

## If some of the high fold changes are getting drowned out due to a large range, you could set a minimum and maximum fold change value
AS_OE_foldchanges <- ifelse(AS_OE_foldchanges > 2, 2, AS_OE_foldchanges)
AS_OE_foldchanges <- ifelse(AS_OE_foldchanges < -2, -2, AS_OE_foldchanges)

cnetplot(ego_AS, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=AS_OE_foldchanges, 
         vertex.label.font=6)

## Repeat for the AC comparision
ids_AC <- bitr(rownames(dge.AC.vp), 
               fromType = "SYMBOL", 
               toType = c("ENSEMBL", "ENTREZID"), 
               OrgDb = "org.Hs.eg.db")

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates_AC <- which(duplicated(ids_AC$SYMBOL) == FALSE)

ids_AC <- ids_AC[non_duplicates_AC, ] 

## Merge the Ensembl IDs with the results     
merged_gene_ids_AC <- merge(x=dge.AC.vp, y=ids_AC, by.x="row.names", by.y="SYMBOL")             

## Extract significant results
AC_sigOE <- subset(merged_gene_ids_AC, padj < 0.05)

AC_sigOE_genes <- as.character(AC_sigOE$ENSEMBL)

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
AC_allOE_genes <- as.character(merged_gene_ids_AC$ENSEMBL)

## Run GO enrichment analysis 
ego_AC <- enrichGO(gene = AC_sigOE_genes, 
                   universe = AC_allOE_genes, 
                   keyType = "ENSEMBL", 
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

## Output results from GO analysis to a table
cluster_summary_AC <- data.frame(ego_AC)

write.csv(cluster_summary_AC, "/Users/ShwetaPipaliya/Documents/PostDoc/projects/viral_transcriptomics/analyses/sars_cov_2/PRJEB43380/enrichment_analyses/clusterProfiler_AC.csv")

## Dotplot gives top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not padj.
dotplot(ego_AC, showCategory=50)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
AC_OE_foldchanges <- AC_sigOE$log2FoldChange

names(AC_OE_foldchanges) <- AC_sigOE$Row.names

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego_AC, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=AC_OE_foldchanges, 
         vertex.label.font=6)

## If some of the high fold changes are getting drowned out due to a large range, you could set a minimum and maximum fold change value
AC_OE_foldchanges <- ifelse(AC_OE_foldchanges > 2, 2, AC_OE_foldchanges)
AC_OE_foldchanges <- ifelse(AC_OE_foldchanges < -2, -2, AC_OE_foldchanges)

cnetplot(ego_AC, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=AC_OE_foldchanges, 
         vertex.label.font=6)

## Repeat for the MS comparision
ids_MS <- bitr(rownames(dge.MS.vp), 
               fromType = "SYMBOL", 
               toType = c("ENSEMBL", "ENTREZID"), 
               OrgDb = "org.Hs.eg.db")

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates_MS <- which(duplicated(ids_MS$SYMBOL) == FALSE)

ids_MS <- ids_MS[non_duplicates_MS, ] 

## Merge the Ensembl IDs with the results     
merged_gene_ids_MS <- merge(x=dge.MS.vp, y=ids_MS, by.x="row.names", by.y="SYMBOL")             

## Extract significant results
MS_sigOE <- subset(merged_gene_ids_MS, padj < 0.05)

MS_sigOE_genes <- as.character(MS_sigOE$ENSEMBL)

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
MS_allOE_genes <- as.character(merged_gene_ids_MS$ENSEMBL)

## Run GO enrichment analysis 
ego_MS <- enrichGO(gene = MS_sigOE_genes, 
                   universe = MS_allOE_genes, 
                   keyType = "ENSEMBL", 
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

## Output results from GO analysis to a table
cluster_summary_MS <- data.frame(ego_MS)

write.csv(cluster_summary_MS, "/Users/ShwetaPipaliya/Documents/PostDoc/projects/viral_transcriptomics/analyses/sars_cov_2/PRJEB43380/enrichment_analyses/clusterProfiler_MS.csv")

## Dotplot gives top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not padj.
dotplot(ego_MS, showCategory=50)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
MS_OE_foldchanges <- MS_sigOE$log2FoldChange

names(MS_OE_foldchanges) <- MS_sigOE$Row.names

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego_MS, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=MS_OE_foldchanges, 
         vertex.label.font=6)

## If some of the high fold changes are getting drowned out due to a large range, you could set a minimum and maximum fold change value
MS_OE_foldchanges <- ifelse(MS_OE_foldchanges > 2, 2, MS_OE_foldchanges)
MS_OE_foldchanges <- ifelse(MS_OE_foldchanges < -2, -2, MS_OE_foldchanges)

cnetplot(ego_MS, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=MS_OE_foldchanges, 
         vertex.label.font=6)

## Repeat for the MC comparision
ids_MC <- bitr(rownames(dge.MC.vp), 
               fromType = "SYMBOL", 
               toType = c("ENSEMBL", "ENTREZID"), 
               OrgDb = "org.Hs.eg.db")

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates_MC <- which(duplicated(ids_MC$SYMBOL) == FALSE)

ids_MC <- ids_MC[non_duplicates_MC, ] 

## Merge the Ensembl IDs with the results     
merged_gene_ids_MC <- merge(x=dge.MC.vp, y=ids_MC, by.x="row.names", by.y="SYMBOL")             

## Extract significant results
MC_sigOE <- subset(merged_gene_ids_MC, padj < 0.05)

MC_sigOE_genes <- as.character(MC_sigOE$ENSEMBL)

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
MC_allOE_genes <- as.character(merged_gene_ids_MC$ENSEMBL)

## Run GO enrichment analysis 
ego_MC <- enrichGO(gene = MC_sigOE_genes, 
                   universe = MC_allOE_genes, 
                   keyType = "ENSEMBL", 
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

## Output results from GO analysis to a table
cluster_summary_MC <- data.frame(ego_MC)

write.csv(cluster_summary_MC, "/Users/ShwetaPipaliya/Documents/PostDoc/projects/viral_transcriptomics/analyses/sars_cov_2/PRJEB43380/enrichment_analyses/clusterProfiler_MC.csv")

## Dotplot gives top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not padj.
dotplot(ego_MC, showCategory=50)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
MC_OE_foldchanges <- MC_sigOE$log2FoldChange

names(MC_OE_foldchanges) <- MC_sigOE$Row.names

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego_MC, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=MC_OE_foldchanges, 
         vertex.label.font=6)

## If some of the high fold changes are getting drowned out due to a large range, you could set a minimum and maximum fold change value
MC_OE_foldchanges <- ifelse(MC_OE_foldchanges > 2, 2, MC_OE_foldchanges)
MC_OE_foldchanges <- ifelse(MC_OE_foldchanges < -2, -2, MC_OE_foldchanges)

cnetplot(ego_MC, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=MC_OE_foldchanges, 
         vertex.label.font=6)

## Repeat for the SC comparison
ids_SC <- bitr(rownames(dge.SC.vp), 
               fromType = "SYMBOL", 
               toType = c("ENSEMBL", "ENTREZID"), 
               OrgDb = "org.Hs.eg.db")

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates_SC <- which(duplicated(ids_SC$SYMBOL) == FALSE)

ids_SC <- ids_SC[non_duplicates_SC, ] 

## Merge the Ensembl IDs with the results     
merged_gene_ids_SC <- merge(x=dge.SC.vp, y=ids_SC, by.x="row.names", by.y="SYMBOL")             

## Extract significant results
SC_sigOE <- subset(merged_gene_ids_SC, padj < 0.05)

SC_sigOE_genes <- as.character(SC_sigOE$ENSEMBL)

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
SC_allOE_genes <- as.character(merged_gene_ids_SC$ENSEMBL)

## Run GO enrichment analysis 
ego_SC <- enrichGO(gene = SC_sigOE_genes, 
                   universe = SC_allOE_genes, 
                   keyType = "ENSEMBL", 
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

## Output results from GO analysis to a table
cluster_summary_SC <- data.frame(ego_SC)

write.csv(cluster_summary_SC, "/Users/ShwetaPipaliya/Documents/PostDoc/projects/viral_transcriptomics/analyses/sars_cov_2/PRJEB43380/enrichment_analyses/clusterProfiler_SC.csv")

## Dotplot gives top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not padj.
dotplot(ego_SC, showCategory=50)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
SC_OE_foldchanges <- SC_sigOE$log2FoldChange

names(SC_OE_foldchanges) <- SC_sigOE$Row.names

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego_SC, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=SC_OE_foldchanges, 
         vertex.label.font=6)

## If some of the high fold changes are getting drowned out due to a large range, you could set a minimum and maximum fold change value
SC_OE_foldchanges <- ifelse(SC_OE_foldchanges > 2, 2, SC_OE_foldchanges)
SC_OE_foldchanges <- ifelse(SC_OE_foldchanges < -2, -2, SC_OE_foldchanges)

cnetplot(ego_SC, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=SC_OE_foldchanges, 
         vertex.label.font=6)

## Perform Gene Set Enrichment Analysis using ClusterProfiler

## Start with the asymptomatic v mild analysis - prepare input
df.am <- read.csv("asymptomatic_v_mild_dge.csv")
colnames(df.am)[1] = "SYMBOL" #rename column 1 as SYMBOL

# remove all values that do not have a padj of 0.05
df.am <- df.am[-which(df.am$padj > 0.05), ]

# convert column 1 into rownames
#df.am.rownames <- df.am %>% remove_rownames %>% column_to_rownames(var="SYMBOL")

# keep only the symbol and the L2FC column
df.am.gsea <- df.am[c("SYMBOL","log2FoldChange")]

# remove row with the empty name 
df.am.gsea <- df.am.gsea[-3614,]

# feature 1: numeric vector
geneList = df.am.gsea[,2]

# feature 2: named vector
names(geneList) = as.character(df.am.gsea[,1])

# omit NA values
geneList_am <- na.omit(geneList_am)

# sort the list in decreasing order (required for clusterProfiler)
sort.am_gene_list = sort(geneList_am, decreasing = TRUE)
#sort.am_gene_list = order(gene_list$log2FoldChange, decreasing = TRUE)

write.csv(sort(geneList_am, decreasing = TRUE), "geneList_GSEAsorted_AM.csv")

# perform gene set enrichment analysis using gseGOvi - This step can take upto 2-3 hours
gse_am <- gseGO(geneList= sort.am_gene_list, 
                ont = "ALL", 
                keyType = "SYMBOL", 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE,
                pAdjustMethod = "none", 
                OrgDb = "org.Hs.eg.db", 
                nPermSimple = 1000)
gse_am
# visualize the results using a dotplot
dotplot(gse_am, showCategory = 15, split=".sign") +
  facet_grid(.~.sign)

# visualize using the encrichment map
emapplot(gse_am, showCategory = 10)

## Repeat for asymptomatic v severe analysis - prepare input
df.as <- read.csv("asymptomatic_v_severe_dge.csv")
colnames(df.as)[1] = "SYMBOL" #rename column 1 as SYMBOL

# remove all values that do not have a padj of 0.05
df.as <- df.as[-which(df.as$padj > 0.05), ]

# convert column 1 into rownames
#df.am.rownames <- df.am %>% remove_rownames %>% column_to_rownames(var="SYMBOL")

# keep only the symbol and the L2FC column
df.as.gsea <- df.as[c("SYMBOL","log2FoldChange")]

# remove row with the empty name 
#df.as.gsea <- df.as.gsea[-3614,]

# feature 1: numeric vector
geneList_as = df.as.gsea[,2]

# feature 2: named vector
names(geneList_as) = as.character(df.as.gsea[,1])

# omit NA values
geneList_as <- na.omit(geneList_as)

# sort the list in decreasing order (required for clusterProfiler)
sort.as_gene_list = sort(geneList_as, decreasing = TRUE)
#sort.as_gene_list = order(gene_list_as$log2FoldChange, decreasing = TRUE)

write.csv(sort(geneList_as, decreasing = TRUE), "geneList_GSEAsorted_AS.csv")

# perform gene set enrichment analysis using gseGOvi
gse_as <- gseGO(geneList= sort.as_gene_list, 
                ont = "ALL", 
                keyType = "SYMBOL", 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE,
                pAdjustMethod = "none", 
                OrgDb = "org.Hs.eg.db", 
                nPermSimple = 1000)
gse_as
# visualize the results using a dotplot
dotplot(gse_as, showCategory = 15, split=".sign") +
  facet_grid(.~.sign)

# visualize using the encrichment map
emapplot(gse_as, showCategory = 10)

## Repeat for asymptomatic v critical analysis - prepare input
df.ac <- read.csv("asymptomatic_v_critical_dge.csv")
colnames(df.ac)[1] = "SYMBOL" #rename column 1 as SYMBOL

# remove all values that do not have a padj of 0.05
df.ac <- df.ac[-which(df.ac$padj > 0.05), ]

# convert column 1 into rownames
#df.ac.rownames <- df.ac %>% remove_rownames %>% column_to_rownames(var="SYMBOL")

# keep only the symbol and the L2FC column
df.ac.gsea <- df.ac[c("SYMBOL","log2FoldChange")]

# remove row with the empty name 
#df.ac.gsea <- df.ac.gsea[-3614,]

# feature 1: numeric vector
geneList_ac = df.ac.gsea[,2]

# feature 2: named vector
names(geneList_ac) = as.character(df.ac.gsea[,1])

# omit NA values
geneList_ac <- na.omit(geneList_ac)

# sort the list in decreasing order (required for clusterProfiler)
sort.ac_gene_list = sort(geneList_ac, decreasing = TRUE)
#sort.ac_gene_list = order(gene_list_ac$log2FoldChange, decreasing = TRUE)

write.csv(sort(geneList_ac, decreasing = TRUE), "geneList_GSEAsorted_AC.csv")

# perform gene set enrichment analysis using gseGOvi
gse_ac <- gseGO(geneList= sort.ac_gene_list, 
                ont = "ALL", 
                keyType = "SYMBOL", 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE,
                pAdjustMethod = "none", 
                OrgDb = "org.Hs.eg.db", 
                nPermSimple = 1000)
gse_ac
# visualize the results using a dotplot
dotplot(gse_ac, showCategory = 15, split=".sign") +
  facet_grid(.~.sign)

# visualize using the encrichment map
emapplot(gse_ac, showCategory = 10)

## Repeat for mild v critical analysis - prepare input
df.mc <- read.csv("mild_v_critical_dge.csv")
colnames(df.mc)[1] = "SYMBOL" #rename column 1 as SYMBOL

# remove all values that do not have a padj of 0.05
df.mc <- df.mc[-which(df.mc$padj > 0.05), ]

# convert column 1 into rownames
#df.mc.rownames <- df.mc %>% remove_rownames %>% column_to_rownames(var="SYMBOL")

# keep only the symbol and the L2FC column
df.mc.gsea <- df.mc[c("SYMBOL","log2FoldChange")]

# remove row with the empty name 
#df.mc.gsea <- df.mc.gsea[-3614,]

# feature 1: numeric vector
geneList_mc = df.mc.gsea[,2]

# feature 2: named vector
names(geneList_mc) = as.character(df.mc.gsea[,1])

# omit NA values
geneList_mc <- na.omit(geneList_mc)

# sort the list in decreasing order (required for clusterProfiler)
sort.mc_gene_list = sort(geneList_mc, decreasing = TRUE)
#sort.mc_gene_list = order(gene_list_mc$log2FoldChange, decreasing = TRUE)

write.csv(sort(geneList_mc, decreasing = TRUE), "geneList_GSEAsorted_MC.csv")

# perform gene set enrichment analysis using gseGOvi
gse_mc <- gseGO(geneList= sort.mc_gene_list, 
                ont = "ALL", 
                keyType = "SYMBOL", 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE,
                pAdjustMethod = "none", 
                OrgDb = "org.Hs.eg.db", 
                nPermSimple = 1000)
gse_mc
# visualize the results using a dotplot
dotplot(gse_mc, showCategory = 15, split=".sign") +
  facet_grid(.~.sign)

# visualize using the encrichment map
emapplot(gse_mc, showCategory = 10)

## Repeat for mild v severe analysis - prepare input
df.ms <- read.csv("mild_v_severe_dge.csv")
colnames(df.ms)[1] = "SYMBOL" #rename column 1 as SYMBOL

# remove all values that do not have a padj of 0.05
df.ms <- df.ms[-which(df.ms$padj > 0.05), ]

# convert column 1 into rownames
#df.ms.rownames <- df.ms %>% remove_rownames %>% column_to_rownames(var="SYMBOL")

# keep only the symbol and the L2FC column
df.ms.gsea <- df.ms[c("SYMBOL","log2FoldChange")]

# remove row with the empty name 
#df.ms.gsea <- df.ms.gsea[-3614,]

# feature 1: numeric vector
geneList_ms = df.ms.gsea[,2]

# feature 2: named vector
names(geneList_ms) = as.character(df.ms.gsea[,1])

# omit NA values
geneList_ms <- na.omit(geneList_ms)

# sort the list in decreasing order (required for clusterProfiler)
sort.ms_gene_list = sort(geneList_ms, decreasing = TRUE)
#sort.ms_gene_list = order(gene_list_ms$log2FoldChange, decreasing = TRUE)

write.csv(sort(geneList_ms, decreasing = TRUE), "geneList_GSEAsorted_MS.csv")

# perform gene set enrichment analysis using gseGOvi
gse_ms <- gseGO(geneList= sort.ms_gene_list, 
                ont = "ALL", 
                keyType = "SYMBOL", 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE,
                pAdjustMethod = "none", 
                OrgDb = "org.Hs.eg.db", 
                nPermSimple = 1000)
gse_ms
# visualize the results using a dotplot
dotplot(gse_ms, showCategory = 15, split=".sign") +
  facet_grid(.~.sign)

## Repeat for severe v critical analysis - prepare input
df.sc <- read.csv("severe_v_critical_dge.csv")
colnames(df.sc)[1] = "SYMBOL" #rename column 1 as SYMBOL

# remove all values that do not have a padj of 0.05
df.sc <- df.sc[-which(df.sc$padj > 0.05), ]

# convert column 1 into rownames
#df.sc.rownames <- df.sc %>% remove_rownames %>% column_to_rownames(var="SYMBOL")

# keep only the symbol and the L2FC column
df.sc.gsea <- df.sc[c("SYMBOL","log2FoldChange")]

# remove row with the empty name 
#df.sc.gsea <- df.sc.gsea[-3614,]

# feature 1: numeric vector
geneList_sc = df.sc.gsea[,2]

# feature 2: named vector
names(geneList_sc) = as.character(df.sc.gsea[,1])

# omit NA values
geneList_sc <- na.omit(geneList_sc)

# sort the list in decreasing order (required for clusterProfiler)
sort.sc_gene_list = sort(geneList_sc, decreasing = TRUE)
#sort.sc_gene_list = order(gene_list_sc$log2FoldChange, decreasing = TRUE)

write.csv(sort(geneList_sc, decreasing = TRUE), "geneList_GSEAsorted_SC.csv")

# perform gene set enrichment analysis using gseGOvi
gse_sc <- gseGO(geneList= sort.sc_gene_list, 
                ont = "ALL", 
                keyType = "SYMBOL", 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE,
                pAdjustMethod = "none", 
                OrgDb = "org.Hs.eg.db", 
                nPermSimple = 1000)
gse_sc
# visualize the results using a dotplot
dotplot(gse_sc, showCategory = 15, split=".sign") +
  facet_grid(.~.sign)

# visualize using the encrichment map
emapplot(gse_sc, showCategory = 10)

## METHOD USING ENTREZIDS
# convert geneIDs to ENTREZIDs
#hs <- org.Hs.eg.db
#my.symbols.am <- c(names(am_gene_list))
#gene_data_am <- AnnotationDbi::select(hs, 
#keys = my.symbols,
#columns = c("ENTREZID", "SYMBOL"), 
#keytype = "SYMBOL")

# merge dataframes
#merge_df.am <- merge(df.am.gsea, gene_data_am, by="SYMBOL")


# TRY USING ENSEMBL
#ids_AM <- bitr(rownames(dge.AM.vp), 
              # fromType = "SYMBOL", 
               #toType = c("ENSEMBL", "ENTREZID"), 
               #OrgDb = "org.Hs.eg.db")

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
#non_duplicates_AM <- which(duplicated(ids_AM$SYMBOL) == FALSE)

#ids_AM <- ids_AM[non_duplicates_AM, ] 

## Merge the Ensembl IDs with the results     
#merged_gene_ids_AM <- merge(x=dge.AM.vp, y=ids_AM, by.x="row.names", by.y="SYMBOL")             

## Extract significant results
#AM_sigOE <- subset(merged_gene_ids_AM, padj < 0.05)

#AM_sigOE_genes <- as.character(AM_sigOE$ENSEMBL)

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
#AM_allOE_genes <- as.character(merged_gene_ids_AM$ENSEMBL)

# extract the log 2 fold change column
#am_gene_list <- AM_sigOE$log2FoldChange

# omit NA values
#am_gene_list <- na.omit(am_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
#sort.am_gene_list = sort(am_gene_list, decreasing = TRUE)

# perform gene set enrichment analysis
#gse_am <- gseGO(geneList= sort.am_gene_list, 
                #ont = "ALL", 
                #keyType = "ENSEMBL", 
                #nPerm = 10000,
                #minGSSize = 3, 
                #maxGSSize = 800, 
                #pvalueCutoff = 0.05, 
                #pAdjustMethod = "BH", 
                #OrgDb = org.Hs.eg.db)

# visualize the results using a dotplot

## Output results from GO analysis to a table
#cluster_summary_AM <- data.frame(ego_AM)

#write.csv(cluster_summary_AM, "/Users/ShwetaPipaliya/Documents/PostDoc/projects/viral_transcriptomics/analyses/sars_cov_2/PRJEB43380/enrichment_analyses/clusterProfiler_AM.csv")


#### CODE BELOW IS NULL = REQUIRES ANNOTATION HUB WHICH I HAVENT FIGURED OUT JUST YET ####
## Visualize enriched terms from the GSE list of DGE genes using DOSE (Disease Ontology Semantic and Enrichment Analysis)
## Explore the grch38 table loaded by the annotables library
#ah = AnnotationHub()
#ah

## Return the IDs for the gene symbols in the DE results
#idx <- grch38$symbol %in% rownames(AM)

#ids <- grch38[idx, ]

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
#non_duplicates <- which(duplicated(ids$symbol) == FALSE)

#ids <- ids[non_duplicates, ] 

## Merge the IDs with the results 
#res_ids <- inner_join(res_tableOE_tb, ids, by=c("gene"="symbol")) 

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
#allOE_genes <- as.character(res_ids$ensgene)

## Extract significant results
#sigOE <- filter(res_ids, padj < 2.9e-06) # use genes that meet the genome wide significance threshold

#sigOE_genes <- as.character(sigOE$ensgene)

### END ###