#!/usr/bin/Rscript --vanilla

## Setup
# install packages from to run subsequent libraries
install.packages(BiocManager)
install.packages(RColorBrewer)
install.packages(pheatmap)
install.packages(ggrepel)
install.packages(devtools)
install.packages(tidyverse)

# install packages from Bioconductor
BiocManager::install("DESeq2")
BiocManager::install("clusterProfiler")
BiocManager::install("DOSE")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("pathview")
BiocManager::install("DEGreport")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("AnnotationHub")
BiocManager::install("ensembldb")

# check library installation
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(clusterProfiler)
library(DEGreport)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(tidyverse)
library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
library(ensembldb)

# print session info
sessionInfo()

## Load in featureCounts data
data <- read.table("data/Mov10_full_counts.txt", header=T, row.names=1) 

meta <- read.table("meta/Mov10_full_meta.txt", header=T, row.names=1)

## Check classes of the data
class(meta)
class(data)

## View data to ensure dataset contains expected samples
View(meta)
View(data)

## plot RNA-Seq featureCount distribution using ggplot2 starting from gene counts with 0 transcripts mapped
ggplot(data) +
   geom_histogram(aes(x = Mov10_oe_1), stat = "bin", bins = 200) + 
   xlim(-5, 500)  +
   xlab("Raw expression counts") +
   ylab("Number of genes")
   
# plot mean versus variance for over-expression replicates
mean_counts <- apply(data[, 3:5], 1, mean)
variance_counts <- apply(data[, 3:5], 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
        geom_point(aes(x=mean_counts, y=variance_counts)) + 
        geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
        scale_y_log10() +
        scale_x_log10()

## use DESeq2 package to normalize count data, library sizes, and RNA composition
vignette("DESeq2")

# Check that sample names match in both files
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

# Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)

# print within script editor
View(counts(dds))

# Perform DESeq2 median of rations method of normalization, size factor, save data matrix
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

## Performing quality control on DE analysis

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA 
plotPCA(rld, intgroup="sampletype")

# Input is a matrix of log transformed values
rld <- rlog(dds, blind=T)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))

# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(meta, pca$x)
ggplot(df) + geom_point(aes(x=PC3, y=PC4, color = sampletype))

# Extract the rlog matrix from the object
rld_mat <- assay(rld)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2

# Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function

head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames

# Plot heatmap
pheatmap(rld_cor)
heat.colors <- brewer.pal(6, "Blues")
pheatmap(rld_cor, color = heat.colors, border_color=NA, fontsize = 10, 
  		fontsize_row = 10, height=20)

# Create DESeq object and run analysis
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
dds <- DESeq(dds)

# Total number of raw and normalized counts per sample
colSums(counts(dds))
colSums(counts(dds, normalized=T))






