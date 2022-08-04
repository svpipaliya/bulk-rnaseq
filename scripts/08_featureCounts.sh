#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH	--cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=shweta.pipaliya@epfl.ch
#SBATCH --job-name=featureCounts_PRJXXXXX
#SBATCH --output=featureCounts_PRJXXXXX_out

mkdir results/counts

module load StdEnv/2020
module load subread/2.0.3

featureCounts -T 6 -s 2 \
  -a /gstore/scratch/hpctrain/chr1_reference_gsnap/chr1_grch38.gtf \
  -o ~/unix_lesson/rnaseq/results/counts/Mov10_featurecounts.txt \
  ~/unix_lesson/rnaseq/results/gsnap/*bam
  
less results/counts/Mov10_featurecounts.txt.summary

less results/counts/Mov10_featurecounts.txt

#end


