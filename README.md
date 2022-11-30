# bulk_transcriptomics
This repo contains scripts and workflows for performing bulk transcriptomics analyses. 

Briefly, the workflow consists of performing quality control steps, read and adapter detection and trimming, read alignment against the human reference genome (build hg38/GRCh38), quality control of the aligned BAM files, sorting and indexing of BAM files, gene-level quantification, differential gene expression analyses, and DEG annotation and protein pathway analyses.

Tools used in the scripts and workflow consists of the following: fastqc, fastp, multiQC, STAR, Salmon, SAMTools, FeatureCounts, DESEQ2, Gene Ontology, STRINGDB, and MSigDB/GSEA.

(in progress) This folder also contains a Snakemake pipeline to perform the abovementioned steps with paired-end read data for automated and reproducible analyses. 
