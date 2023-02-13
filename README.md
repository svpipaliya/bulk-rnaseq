# Bulk Transcriptomics for Human RNASeq Data
This repo contains bash scripts and a snakemake workflow for performing bulk transcriptomic analyses, specifically using samples originating from human tissue/cell culture. Workflow can be used for sequencing reads in paired or single-end form generated using Illumina HiSeq/MiSeq or BGI Genomics.  

 - Briefly, the steps consist of  quality control to assess read contamination and base quality, adapter detection and trimming (Illumina only; please provide adapter sequences in-case using another technology), read alignment against the human reference genome (build hg38/GRCh38), quality control of the aligned BAM files, sorting and indexing of BAM files, gene-level quantification, differential gene expression analyses, and DEG annotation and protein pathway analyses.

Tools used:
 -   fastqc (QC)
 -   fastp (Adapter removal and trimming)
 -   multiQC (Data visualization for QC'd results)
 -   STAR (Genome-alignment to BAM format)
 -   Salmon (Transcriptome-alignment and transcript-level quantification)
 -   SAMTools (BAM sorting and indexing)
 -   FeatureCounts (Gene-level quantification)
 -   tximport (Conversion of Salmon transcript-level abundances to gene-level counts)
 -   DESEQ2 (Count normalization and differential gene expression analyses)
 -   Clusterprofiler/GSEA/Gene Ontology (gene set enrichment analyses)
 -   STRINGDB (protein network reconstruction)

The Snakemake workflow is still under progress and contents will be modified frequently. Below is a DAG of the rules used in the workflow.

![Screenshot 2023-01-23 at 17 31 11](https://user-images.githubusercontent.com/61172011/214095023-591e9fc1-dff0-4798-ac86-416f29dfc44c.png)

An additional README section will be added to run the snakemake workflow in an slurm HPC environment as well as instructions on set up using Conda/Mamba.

Author: Shweta Pipaliya

Date Updated: 13.2.2023
