#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH	--cpus-per-task=8
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-user=shweta.pipaliya@epfl.ch
#SBATCH --job-name=samtools_sort
#SBATCH --output=samtools_sort.o
#SBATCH --error=samtools_sort.e

## load modules
module load gcc/8.4.0
module load samtools/1.10

## Set directory for bam files
SET DIR = /work/backup/gr-fe/yassine/COVID19_GSE171110/GSE171110/02_aligned_reads/

## sort star output bam alignments by read name using samtools
for f in $DIR/*.bam; do  samtools sort -n -@ 7 $f -o ${f/.bam/name_sorted.bam}; done

## end
