#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH	--cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=shweta.pipaliya@epfl.ch
#SBATCH --job-name=fastq_dump_PRJNA2722046
#SBATCH --output=fastq_dump_PRJNA2722046_out

module load gcc/9.3.0
module load sra-toolkit/2.10.8
module load nixpkgs/16.09
module load fastqc/0.11.9
module load star/2.7.9a
module load StdEnv/2020
module load subread/2.0.3
module spider r/4.0.2

# for paired end reads only
fastq-dump --split-3  $1

# for every SRR in the list of SRRs file
for srr in $(cat SRR_Acc_List.txt)
  do
  # call the bash script that does the fastq dump, passing it the SRR number next in file
  sbatch 01_inner_fastq_dump.slurm $srr
  sleep 1	# wait 1 second between each job submission

#  Running FASTQC
fastqc -t 6 *.fq


# make a results directory
mkdir /results/fastqc

# Moving files to our results directory
mv *fastqc* ../results/fastqc/

#parameter specification for running STAR *** replace the directory paths ***
STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir hg38_index \
--genomeFastaFiles Homo_sapiens.GRCh38.fa \
--sjdbGTFfile /Homo_sapiens.GRCh38.92.gtf \
--sjdbOverhang 99

mkdir results/counts

module load StdEnv/2020
module load subread/2.0.3

featureCounts -T 6 -s 2 \
  -a /gstore/scratch/hpctrain/chr1_reference_gsnap/chr1_grch38.gtf \
  -o ~/unix_lesson/rnaseq/results/counts/Mov10_featurecounts.txt \
  ~/unix_lesson/rnaseq/results/gsnap/*bam
  
less results/counts/Mov10_featurecounts.txt.summary

less results/counts/Mov10_featurecounts.txt

# run DESeq2 R-Script
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
   
   mean_counts <- apply(data[, 3:5], 1, mean)
variance_counts <- apply(data[, 3:5], 1, var)
df <- data.frame(mean_counts, variance_counts)

# plot mean versus variance for over-expression replicates
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

