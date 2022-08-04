#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH	--cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=shweta.pipaliya@epfl.ch
#SBATCH --job-name=STAR_genome_index
#SBATCH --output=STAR_genome_index_out

# load all modules for runnning the STAR v.2.7.9a aligner
module load StdEnv/2020
module load star/2.7.9a

# parameter specification for running STAR to generate a genome index
STAR --runThreadN 16 \
--runMode genomeGenerate \
--genomeDir /project/6001383/pip17/sars_cov_2/data/GSE172114/hg38_ref/ \
--genomeFastaFiles /project/6001383/pip17/sars_cov_2/data/GSE172114/hg38_ref/GCA_000001405.15_GRCh38_full_analysis_set.fasta \
--sjdbGTFfile /project/6001383/pip17/sars_cov_2/data/GSE172114/hg38_ref/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf \
--sjdbOverhang 99

#end