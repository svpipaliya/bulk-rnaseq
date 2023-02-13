# Bulk Transcriptomics for Human RNASeq Data
This repo contains bash scripts and a snakemake workflow for performing bulk transcriptomic analyses, specifically using samples originating from human tissue/cell culture. Workflow can be used for sequencing reads in paired or single-end form generated using Illumina HiSeq/MiSeq or BGI Genomics.  

 - Briefly, the steps consist of  quality control to assess read contamination and base quality, adapter detection and trimming (Illumina only; please provide adapter sequences in-case using another technology), read alignment against the human reference genome/transcriptome (build hg38/GRCh38), quality control of the aligned BAM files, sorting and indexing of BAM files, gene-level quantification, differential gene expression analyses, and DEG annotation and protein pathway analyses.

Tools and dependencies:
 -   pandas = 0.23
 -   star=2.7.0
 -   fastp=0.20.1-0
 -   subread=2.0.1-0
 -   multiqc=1.9-0
 -   yaml = 0.2.5
 -   fastqc =0.11.9=0
 -   samtools=1.3.1
 -   salmon=1.4.0
 -   tximport
 -   DESEQ2
 -   Clusterprofiler
 -   GSEA
 -   Gene Ontology
 -   STRINGDB

The Snakemake workflow is still under progress and contents will be modified frequently. Below is the current DAG of the rules used in the workflow.

![Screenshot 2023-01-23 at 17 31 11](https://user-images.githubusercontent.com/61172011/214095023-591e9fc1-dff0-4798-ac86-416f29dfc44c.png)

An additional README section will be added to run the snakemake workflow in an slurm HPC environment as well as instructions on set up using Conda/Mamba.

Author: Shweta Pipaliya

Date Updated: 13.2.2023
