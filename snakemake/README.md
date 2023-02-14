# Description
This page contains instructions for running this snakemake workflow on your slurm HPC. The pipeline is used to process single or paired-end reads and generates a raw gene and transcript-level count matrix that can be used for differential gene expression analysis (implementation of differential gene expression and enrichment analysis steps are currently underway). Additional outputs include trimmed fastq files, `.bam` . 
Currently the pipeline is best-suited for Illumina sequencing data. However, read-processing steps can be modified by the user by providing own adapter sequences. 

## Pipeline overview
1. Raw fastq files (.gz compressed or uncompressed) will be checked for quality using `fastqc`
2. Trimming and adapter removal will be peformed using `fastp`
3. 
