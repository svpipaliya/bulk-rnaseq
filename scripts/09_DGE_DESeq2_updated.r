### (c) Shweta V. Pipaliya 
### 2022-15-05 
### This script is used to perform and visualize differential gene expression analyses of the gene count data using DESeq2
### Input files required: count matrix and sample metadata file
### How to run this script from command-line: Rscript --vanilla 3.dge.svp.R
### Bioconductor packages in this script require R v.4.2.0

## Set seed 
set.seed(777)

## Set current working directory to output DGE summary results and plots
wd <- getwd() 
if (!is.null(wd)) setwd(wd) # set current directory as working directory; wd should have the two input files and this rscript
dir.create("dge.results.svp") # create a new results folder 
list.files(wd)

## Install and load required libraries
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

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

### FIRST PART OF THE SCRIPT: Data formatting, sample and gene qc, count normalization, and DGE using DESeq2 ###

## Step 1: load gene counts for all genes from healthy and disease samples as a matrix
counts <- as.matrix(read.csv("GSE171110_counts_rs.csv", header = TRUE, row.names = "Geneid", sep=","))
counts[is.na(counts)] <- 0 # convert any NA values to 0 in the counts matrix

## Step 2: load sample metadata with individual samples IDs corresponding to healthy vs. disease 
metaData <- read.csv("metadata.deseq.csv", row.names = 1, sep=",")

## Step 3: ensure all value names in the count data match the order of column names in the metaData csv
all(colnames(counts) %in% rownames(metaData))
all(colnames(counts) == rownames(metaData))

#which(colnames(counts) != rownames(metaData)) # check which rownames and colnames are not right if returned as false
#colnames(counts)[46]
#rownames(metaData)[46]

## Step 4: create object of class DESeqDataSet to create a matrix that can be read by DESeq2
deseq2Data <- DESeqDataSetFromMatrix(countData = counts, 
                                     colData = metaData,
                                     design = ~ sample_group)
print(deseq2Data)

## Step 5: transform gene counts for PCA for sample-level QC 
pdf("dge.results.svp/plot.GSE172114.sample.pca.pdf")

rld <- rlog(deseq2Data, blind=TRUE)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat)) # run prcomp for PCA
df.pca = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], 
                    PC3 = pca$x[,3], PC4 = pca$x[,4]) 

plot.pca <- plotPCA(rld,intgroup="sample_group") # use DESeq2's plotPCA() for plotting
plot.pca + ggtitle("Principal Component Analysis: Visualizing sample variance")

dev.off()

## Step 6: remove disease sample outliers clustering away (PC1 > 10) 
#outliers <- as.character(subset(colnames(deseq2Data), df.pca$PC1 > 10))
#print (outliers)

#deseq2Data.sqc <- deseq2Data[, !(colnames(deseq2Data) %in% outliers)] 

## Step 7: extract the rlog matrix from the QC'd deseq2 data object
rld.sqc <- rlog(deseq2Data, blind=T)
rld_mat.sqc <- assay(rld.sqc)   

## Step 8: compute pairwise correlation values for post-QC samples and plot a heatmap
pdf("dge.results.svp/plot.sampleQC.heatmap.pdf")

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

## Step 9: perform gene-level QC by filtering out rows with less than 5 reads across all the samples
dim(deseq2Data)
dim(deseq2Data[rowSums(counts(deseq2Data)) > 5, ]) # assess filter effect
deseq2Data.gqc <- deseq2Data[rowSums(counts(deseq2Data)) > 5, ] # run and output filtered data

## Step 10: set 'Healthy' as reference for comparison against 'Disease' using relevel
deseq2Data.gqc$sample_group <- relevel(deseq2Data.gqc$sample_group, ref = "healthy")
print(deseq2Data.gqc$sample_group) # important: levels should be healthy first then disease

## Step 11: perform DGE analyses on the filtered dataset using the DESeq() on multiple cores - may take several minutes
deseq2Data_DGE <- DESeq(deseq2Data.gqc, parallel = TRUE)
deseq2Data_r0.05 <- results(deseq2Data_DGE, alpha = 0.05)
dgeResults <- results(deseq2Data_DGE)

head(dgeResults)

dgeSum <- summary(dgeResults) # view summary summary statistics for up- and down-regulated genes
dgeSum <- summary(deseq2Data_r0.05)

## Step 12: sort and write out DGE results for all genes by adjusted p-values (Benjamini-Hochberg FDR method)
dgeResults <- dgeResults[order(dgeResults$padj),]
head(dgeResults)

deseq2Data_r0.05 <- deseq2Data_r0.05[order(deseq2Data_r0.05$padj),]

write.csv(as.data.frame(deseq2Data_r0.05[order(deseq2Data_r0.05$padj),]), file = "dge.results.svp/all.nc_v_critical_dge.csv")
#write.csv(as.data.frame(dgeResults[order(dgeResults$padj),] ), file="dge.results.svp/all.nc_v_critical_dge.csv")

## Step 13: get gene names for the top 50 dge
top50_sigOE_genes <- rownames(dgeResults[1:50, ]) # get rownames
top50_sigOE_dge <- dgeResults[1:50, ] # extract dge results for the top 50 genes

write.csv(as.data.frame(dgeResults[1:50, ] ), file="dge.results.svp/top50.nc_v_critical_dge.csv")

## Step 13: get and write out normalized counts for all genes
normalized_counts <- counts(deseq2Data_DGE, normalized=TRUE)
head(normalized_counts)

write.csv(as.data.frame(counts(deseq2Data_DGE, normalized=TRUE) ), file="dge.results.svp/all.nc_v_critical_nc.csv")

## Step 14: get and write out normalized counts for the top 50 differentially expressed genes
top50_sigOE_norm <- normalized_counts[top50_sigOE_genes, ]

write.csv(as.data.frame(normalized_counts[top50_sigOE_genes, ] ), file="dge.results.svp/top50.nc_v_critical_nc.csv")

### SECOND PART OF THE SCRIPT: Visualize log fold change and normalized counts for biomarkers of interest ###

## Step 15: generate an MA plot for all p-adj DGE results with biomarkers highlighted
pdf(file = "dge.results.svp/plot.MA.pdf") 

df.all.deSeqRes <- data.frame(dgeResults, stringsAsFactors = FALSE)
df.all.deSeqRes <- tibble::rownames_to_column(df.all.deSeqRes, "row_names") # convert original dgeResults df with gene row names as a separate column for MA plotting
colnames(df.all.deSeqRes)[1] <- "Gene" # change column 1name to Gene

df.all.deSeqRes$Significance <- ifelse(df.all.deSeqRes$padj <= 0.05, "Yes", "No") # specify if gene significance meets padj =< 0.05 threshold
df.all.deSeqRes$genelabels <- factor(df.all.deSeqRes$Gene, levels = c(top50_sigOE_genes)) # add label column for biomarkers

ggplot(df.all.deSeqRes) +
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
  ggtitle("MA plot: Log fold-change vs. mean expression") +
  theme(plot.title = element_text(size = rel(1.0), hjust = 0.5),
        axis.title = element_text(size = rel(0.75))) 

dev.off()

## Step 16: generate a volcano plot depicting fold change (>=1.5) of all genes with biomarkers highlighted
pdf(file = "dge.results.svp/plot.volcano.pdf") 

fc_all.deSeqRes <- df.all.deSeqRes  %>%
  mutate(Fold.change_1.5 = padj <= 0.05 & abs(log2FoldChange) >= 0.58)
fc_all.deSeqRes$genelabels <- factor(df.all.deSeqRes$Gene, levels = c(top50_sigOE_genes)) # add label column for biomarkers

ggplot(fc_all.deSeqRes, aes(x = log2FoldChange, 
                            y = -log10(padj))) +
  geom_point(aes(colour = Fold.change_1.5), 
             size = 0.25) +
  geom_label_repel(aes(label = genelabels),
                   box.padding = unit(0.35, "lines"),
                   point.padding = 0.2,
                   nudge_x = 0.0, 
                   nudge_y = 0.5, 
                   size = 2) + 
  ggtitle("Volcano plot: depiction of gene expression levels") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  xlim(-7, 7) +
  ylim(0, 10) +
  theme(legend.title = element_text(size = rel(0.25)),
        legend.text = element_text(size = rel(0.25)),
        plot.title = element_text(size = rel(1.00), 
                                  hjust = 0.5),
        axis.title = element_text(size = rel(0.25))) 

dev.off()

## Step 17: plot boxplots of normalized count comparisons for top 50 biomarkers in healthy vs. disease samples
pdf("dge.results.svp/plot.boxplot.pdf")

goi_NC <- read.csv("dge.results.svp/top50.nc_v_critical_nc.csv") # read in the normalized counts csv for biomarkers from previous steps
metaData_plot <- read.csv("metadata.deseq.csv") # read in the metadata for table joining 

colnames(goi_NC)[1] <- "gene" # change column 1 name to gene for plotting
gathered_goi_NC <- goi_NC %>% # gather normalized counts for each sample (cols 2-6) into a single column
  gather(colnames(goi_NC) [2:55], key = "sample", value = "normalized_counts")

gathered_goi_NC <- inner_join(metaData_plot, gathered_goi_NC)  # join metadata (healthy vs. disease) with the gathered df

ggplot(gathered_goi_NC, 
       aes(x = gene, 
           y = normalized_counts)) +
  scale_y_log10() +
  ylab("Normalized Counts") +
  xlab("Gene") +
  ggtitle("Box plot: Normalized counts for biomarkers of interest") + 
  geom_boxplot(aes(fill = sample_group), 
               position = position_dodge(0.9)) +
  facet_wrap(~ gene, scales = "free")

dev.off()

## Step 18 (last): draw a heatmap for Z-scores of the log-scaled normalized gene counts for the biomarkers
pdf(file="dge.results.svp/plot.nc.heatmap.pdf")

heatmap_data <- goi_NC %>%
  select(1:55) %>%
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