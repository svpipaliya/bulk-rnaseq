#!/bin/csh
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH	--cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-user=shweta.pipaliya@epfl.ch
#SBATCH --job-name=featureCounts_GSE171110
#SBATCH --output=featureCounts_GSE171110.o
#SBATCH --error=featureCounts_GSE171110.e

## load modules
module load gcc/8.4.0
module load subread/2.0.0

## set variables for reference datasets and directory
set REF_GTF = /work/backup/gr-fe/pipaliya/hg38/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf
set BAM_DIR = /work/backup/gr-fe/pipaliya/transcriptomics/datasets/sars_cov2/GSE171110/04_samtools_sort/*.bam 

## make output directory
#mkdir fc_GSE171110.output

## run FeatureCounts on the BAM files to generate gene counts for reads
featureCounts -a $REF_GTF \
 -s 2 \
 -t exon \
 -T 8 \
 -g gene_id \
 -o GSE171110_counts_rs.txt \
 $BAM_DIR
  
## end



