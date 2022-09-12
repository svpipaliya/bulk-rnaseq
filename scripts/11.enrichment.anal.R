## This script is used to visualize enriched terms from the GSE list of DGE genes using DOSE (Disease Ontology Semantic and Enrichment Analysis)
### (c) Shweta V. Pipaliya 
### 2022-12-09 
### Input files required:geneList
### How to run this script from command-line: Rscript --vanilla 11.enrichment.anal.R
### Bioconductor packages in this script require R v.4.2.0

## Set seed 
set.seed(777)

## Set current working directory to output DGE summary results and plots
wd <- getwd() 
if (!is.null(wd)) setwd(wd) # set current directory as working directory; wd should have the two input files and this rscript
dir.create("enrichment.results") # create a new results folder 
list.files(wd)

## Install and load required libraries
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)

#if (!require(ggplot2)) install.packages('ggplot2')
library(ggplot2)

#BiocManager::install("DOSE")
library(DOSE)

#BiocManager::install("biomaRt")
library(biomaRt)

#BiocManager::install("enrichplot")
library(enrichplot)

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

install.packages("ggnewscale")
library(ggnewscale)

install.packages("ggupset")
library(ggupset)

# map the gene names to get entrezgene_ids

mapping <- getBM(
  attributes = c('entrezgene_id', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values = geneList_AM,
  mart = hsmart
)

hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl") # first convert the gene IDs for ENSEMBL and ENTREZ ID
#mapping_char <- as.character(mapping$entrezgene_id) #convert the values to character instead of numeric for the DOSE package
write.csv(as.data.frame(mapping), file="/Users/ShwetaPipaliya/Documents/PostDoc/projects/viral_transcriptomics/analyses/sars_cov_2/GSE172114/06_DGE/dge.results.svp/entrez_id.csv")

# read in the geneList data
gse <- read.csv("/Users/ShwetaPipaliya/Documents/PostDoc/projects/viral_transcriptomics/analyses/sars_cov_2/GSE172114/07_enrichment_pathway_analyses/geneList.csv", header = TRUE, sep=",") # use the normalized counts csv for the GSE set of genes
# feature 1: numeric vector
geneList = gse[,2]
## feature 2: named vector
names(geneList) = as.character(gse[,1])
## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)

de <- names(geneList)[abs(geneList) > 2]

edo <- enrichDGN(de)

# use the DOSE and the enrichPlot package + entrez gene IDs mapped to visualize the results
edo2 <- gseDO(geneList, 
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH",
              by = "fgsea")
dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")

# perform gene-network concept analysis
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))


p1 <- cnetplot(edox, node_label="category", 
               cex_label_category = 1.2) 
p2 <- cnetplot(edox, node_label="gene", 
               cex_label_gene = 0.8) 
p3 <- cnetplot(edox, node_label="all") 
p4 <- cnetplot(edox, node_label="none", 
               color_category='firebrick', 
               color_gene='steelblue') 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

#heatmap like functional classification
p1 <- heatplot(edox, showCategory=5)
p2 <- heatplot(edox, foldChange=geneList, showCategory=5)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])


# generate an upset plot
upsetplot(edo)

### END ###