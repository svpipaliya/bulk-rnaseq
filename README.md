# bulk_transcriptomics
This repo contains scripts and workflows for performing bulk transcriptomics analyses, specifically using samples originating from human tissue/cell culture. Workflow can be used for sequencing reads (paired or single-end) from Illumina HiSeq/MiSeq or BGI Genomics.  

Briefly, the steps consist of  quality control to assess read contamination and base quality, adapter detection and trimming, read alignment against the human reference genome (build hg38/GRCh38), quality control of the aligned BAM files, sorting and indexing of BAM files, gene-level quantification, differential gene expression analyses, and DEG annotation and protein pathway analyses.

Tools used are as following: fastqc (Quality Control), fastp (read-trimming), multiQC (Data visualization for QC'd results), STAR (alignment to BAM format), Salmon (alignment and transcript-level quantification), SAMTools (sorting and indexing of BAM format files), FeatureCounts (in case Salmon is not used for gene-level quantification), DESEQ2 (normalization of counts and differential gene expression analyses), MSigDB/GSEA/Gene Ontology (gene set enrichmenet analyses), and STRINGDB (protein network reconstruction).


(in progress) This folder also contains a Snakemake pipeline to perform the abovementioned steps with paired-end read data. Below is a DAG of the rules used in the workflow.

![Screenshot 2023-01-23 at 17 31 11](https://user-images.githubusercontent.com/61172011/214095023-591e9fc1-dff0-4798-ac86-416f29dfc44c.png)


Author: Shweta Pipaliya

Date Updated: 23.1.2023
