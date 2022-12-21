# bulk_transcriptomics
This repo contains scripts and workflows for performing bulk transcriptomics analyses, specifically using samples originating from human tissue/cell culture. 

Briefly, the workflow consists of performing quality control steps, read and adapter detection and trimming, read alignment against the human reference genome (build hg38/GRCh38), quality control of the aligned BAM files, sorting and indexing of BAM files, gene-level quantification, differential gene expression analyses, and DEG annotation and protein pathway analyses.

Tools used in the scripts and workflow consists of the following: fastqc (Quality Control), fastp (read-trimming), multiQC (Data visualization for QC'd results), STAR (alignment to BAM format), Salmon (alignment and transcript-level quantification), SAMTools (sorting and indexing of BAM files), FeatureCounts (gene-level quantification), DESEQ2 (differential gene expression analyses), MSigDB/GSEA/Gene Ontology (gene set enrichmenet analyses), and STRINGDB (protein network reconstruction).

(in progress) This folder also contains a containerized Snakemake pipeline to perform the abovementioned steps with paired-end read data in an automated analyses.

Author: Shweta Pipaliya

Date Updated: 30.11.2022
