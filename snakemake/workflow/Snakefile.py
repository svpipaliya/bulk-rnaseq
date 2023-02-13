## Author: Shweta Pipaliya ##
## Date: 09.01.2023 ##
## Snakefile to run the RNASEQ workflow ##
## for a dry-run use: snakemake --dryrun ##
## to generate a DAG of the steps use: snakemake --dag ##

##### libraries for data processing #####
import os
import glob
import pandas as pd

##### read in input consisting of accession and fw/rv read wildcards #####
fibro,FRR = glob_wildcards("../input/rawReads/{fibro}_{frr}.fastq.gz")

##### Load RNASeq Rules #####
include: "rules/1.index.smk"
include: "rules/2.trim.smk"
include: "rules/3.align.smk"
include: "rules/4.alignqc.smk"
include: "rules/5.counts.smk"

#### Target Rules #####
rule all: 
	input: 
		expand("rawQC/{fibro}_{frr}_fastqc.{extension}", fibro=fibro, frr=FRR, extension=["gz","html"]),
		directory("../resources/salmon_hg38_index"),
		directory("../resources/star_genome"),
		"../output/bamQC/multiqc_star_report.html",
		"../output/featureCounts/hgid_feature_counts.txt", 
		#expand("../output/salmonQuant/{fibro}_quant.sf", fibro=fibro)

### End - Run normalisation and DGE steps using the DESEQ2 Rscript seperately ###
