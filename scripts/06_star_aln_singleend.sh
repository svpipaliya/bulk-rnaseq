#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH	--cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=shweta.pipaliya@epfl.ch
#SBATCH --job-name=STAR_alignment_SRR14240795
#SBATCH --output=STAR_alignment_SRR14240795_out

# load all modules for runnning the STAR v.2.7 aligner
module load gcc
module load star/2.7.9a

# parameter specification for running STAR

STAR --genomeDir /home/pip17/scratch/hg38/ --runThreadN 16 --readFilesIn /home/pip17/scratch/sars_cov_2/data/GSE172114/filtered_pe_reads/filtered_reads_fastp_default/SRR14240795_trim_1.fastq.gz /home/pip17/scratch/sars_cov_2/data/GSE172114/filtered_pe_reads/filtered_reads_fastp_default/SRR14240795_trim_2.fastq.gz --sjdbGTFfile /home/pip17/scratch/hg38/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --quantMode GeneCounts --outFileNamePrefix SRR14240795_map

# end
