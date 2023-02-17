# Bulk Transcriptomics for Human RNASeq Data
This repo contains bash scripts and a snakemake workflow for performing bulk transcriptomic analyses, specifically using samples originating from human tissue/cell culture. Workflow can be used for sequencing reads in paired or single-end form generated using Illumina HiSeq/MiSeq or BGI Genomics.  

 - Briefly, the steps consist of  quality control to assess read contamination and base quality, adapter detection and trimming (Illumina only; please provide adapter sequences in-case using another technology), read alignment against the human reference genome/transcriptome (build hg38/GRCh38), quality control of the aligned BAM files, sorting and indexing of BAM files, gene-level quantification, differential gene expression analyses, and DEG annotation and protein pathway analyses.

Tools and dependencies:
 -   [snakemake=7.21.0](https://snakemake.readthedocs.io/en/v7.21.0/)
 -   [snakemake-minimal=7.21.0](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
 -   [python=3.11.0](https://docs.python.org/3/whatsnew/3.11.html)
 -   [pandas = 0.23](https://pandas.pydata.org/)
 -   [star=2.7.0](https://github.com/alexdobin/STAR)
 -   [fastp=0.20.1-0](https://github.com/OpenGene/fastp)
 -   [subread=2.0.1-0](https://subread.sourceforge.net/)
 -   [multiqc=1.9-0](https://multiqc.info/)
 -   [yaml = 0.2.5](https://yaml.org/)
 -   [fastqc =0.11.9=0](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
 -   [samtools=1.3.1](http://www.htslib.org/)
 -   [salmon=1.4.0](https://salmon.readthedocs.io/en/latest/salmon.html)
 -   [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html)
 -   [DESEQ2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
 -   [Clusterprofiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
 -   [GSEA](https://www.gsea-msigdb.org/gsea/index.jsp)
 -   [Gene Ontology](http://geneontology.org/)
 -   [STRINGDB](https://string-db.org/)

The Snakemake workflow is still under progress and contents will be modified frequently. Below is the current DAG of the rules used in the workflow.

![Screenshot 2023-01-23 at 17 31 11](https://user-images.githubusercontent.com/61172011/214095023-591e9fc1-dff0-4798-ac86-416f29dfc44c.png)

An additional [README](https://github.com/svpipaliya/bulk-rnaseq/tree/main/snakemake) section is available to run the snakemake workflow in an slurm HPC environment as well as instructions on set up using Conda/Mamba.

Author: Shweta Pipaliya

Date Updated: 17.2.2023
