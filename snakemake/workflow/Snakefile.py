## Author: Shweta Pipaliya ##
## Date: 09.01.2023 ##
## Snakefile to run the RNASEQ workflow ##
## for a dry-run use: snakemake --dryrun ##
## to generate a DAG of the steps use: snakemake --dag ##

##### Libraries for data processing #####
import os
import glob
import pandas as pd
import sys

##### Read in input consisting of accession and fw/rv read wildcards #####
fibro,FRR = glob_wildcards("../input/rawReads/{fibro}_{frr}.fastq.gz")

##### Add code here to read .tsv/csv containing sample info + metadata using pd #####
units = pd.read_table(



##### Create output directories ######
try:
	dirs = ['trimmedReads', 'starAligned', 'multiQC', 'salmonIndex',
			'salmonAligned', 'samsort', 'featureCounts', 'salmonQuant']		
	for items in dirs:
		os.mkdir(items)
except FileExistsError:
	pass

##### Target Rules ######
rule all: 
	input: 
		expand("rawQC/{fibro}_{frr}_fastqc.{extension}", fibro=fibro, frr=FRR, extension=["zip","html"]),
		directory("../resources/star2.7_hg38_index"),
		#directory("../resources/salmon_hg38_index"),
		directory("../resources/star_genome"),
		"../output/bamQC/multiqc_star_report.html",
		"../output/featureCounts/hgid_feature_counts.txt", 
		#expand("../output/salmonQuant/{fibro}_quant.sf", fibro=fibro)
		
##### Load RNASeq Rules #####
include: "rules/starIndex.smk"
#include: "rules/salmonIndex.smk"
include: "rules/fastQC.smk"	
include: "rules/trim.smk"
include: "rules/starAlign.smk"
#include: "rules/salmonAlign.smk"	
include: "rules/alignqc.smk"
include: "rules/bamsort.smk"
include: "rules/featureCounts.smk"
#include: "rules/salmonQuant.smk"

##### End of Snakemake Run - Perform count normalisation and differential gene expression steps using the DESEQ2 Rscript seperately ######
