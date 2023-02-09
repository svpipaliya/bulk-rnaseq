### (c) Shweta V. Pipaliya 
### 2023-02-09 
## This script converts Salmon quant.sf into gene-level counts 
### run from command-line: Rscript --vanilla salmon_tximport.svp.R
### Bioconductor packages in this script require R v.4.2.0

## Set seed 
set.seed(777)

## Set current working directory to output DGE summary results and plots
wd <- getwd() 
if (!is.null(wd)) setwd(wd) # set current directory as working directory; wd should have the two input files and this rscript
dir.create("tximport_salmon.results.svp") # create a new results folder 
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

#BiocManager::install("BiocParallel")
library(BiocParallel)
register(MulticoreParam(4)) # register the number of nodes on multicoreParam for DESeq2 run

#if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)

install.packages("devtools")
devtools::install_github("stephenturner/annotables")

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
txi <- tximport(samples, type="salmon", 
                tx2gene=tx2gene[,c("refseq_mrna", "hgnc_symbol")], 
                countsFromAbundance="lengthScaledTPM", ignoreTxVersion=TRUE)
attributes(txi) # read the output from the txi object 

head(txi$counts)

## end
