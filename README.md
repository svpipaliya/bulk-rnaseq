# Bulk Transcriptomics for Human RNASeq Data
This repo contains scripts and workflows for performing bulk transcriptomics analyses, specifically using samples originating from human tissue/cell culture. Workflow can be used for sequencing reads (paired or single-end) from Illumina HiSeq/MiSeq or BGI Genomics.  

 - Briefly, the steps consist of  quality control to assess read contamination and base quality, adapter detection and trimming, read alignment against the human reference genome (build hg38/GRCh38), quality control of the aligned BAM files, sorting and indexing of BAM files, gene-level quantification, differential gene expression analyses, and DEG annotation and protein pathway analyses.
- Tools used
-  fastqc (Quality Control), fastp (read-trimming), multiQC (Data visualization for QC'd results), STAR (alignment to BAM format), Salmon (alignment and transcript-level quantification), SAMTools (sorting and indexing of BAM format files), FeatureCounts (alternative to Salmon for gene-level quantification), tximport (conversion of transcript-level abundances to gene-level counts), DESEQ2 (normalization of counts and differential gene expression analyses), MSigDB/GSEA/Gene Ontology (gene set enrichmenet analyses), and STRINGDB (protein network reconstruction).


- This repository consists of a Snakemake pipeline to perform the abovementioned steps with paired-end read data. Below is a DAG of the rules used in the workflow.

![Screenshot 2023-01-23 at 17 31 11](https://user-images.githubusercontent.com/61172011/214095023-591e9fc1-dff0-4798-ac86-416f29dfc44c.png)

- Currently, updates are still in progress and contents will be modified frequently. An additional README section will be added to run the snakemake workflow in an HPC environment and installation using Conda.

Author: Shweta Pipaliya

Date Updated: 23.2.2023
