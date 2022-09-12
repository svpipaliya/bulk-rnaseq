#!/bin/bash
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH	--cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=shweta.pipaliya@epfl.ch
#SBATCH --job-name=STAR_aln_PRJEB43380
#SBATCH --output=STAR_aln_PRJEB43380_out

# load all modules for runnning the STAR v.2.7 aligner
module load StdEnv/2020
module load star/2.7.9a

# Define variables
GNM_DIR=/project/6001383/pip17/sars_cov_2/data/hg38_ref

# Path to the trimmed reads
FQ_DIR=/project/6001383/pip17/sars_cov_2/data/PRJEB43380/02_trimmed_reads

# Directory for output files
#BAM_DIR=/project/6001383/pip17/sars_cov_2/data/PRJEB43380/03_star_alignment

# for every filename in the list, run STAR alignment steps using the following parameters
for base in $(cat Acc_list_PRJEB43380.txt)
  do
  echo -e $base
  
  # define R1 and R2 fastq filenames
  fq1=$FQ_DIR/${base}_trim_1.fq.gz
  fq2=$FQ_DIR/${base}_trim_2.fq.gz
 
     STAR --genomeDir $GNM_DIR \
     --runThreadN 1 \
     --readFilesIn $fq1 $fq2 \
     --readFilesCommand zcat \
     --sjdbGTFfile $GNM_DIR/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf \
     --outSAMtype BAM Unsorted SortedByCoordinate \
     --outSAMunmapped Within \
     --outSAMattributes Standard \
     --outFileNamePrefix $base"_"
done
​
# end