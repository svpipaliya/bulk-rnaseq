#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH	--cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=shweta.pipaliya@epfl.ch
#SBATCH --job-name=fastq_dump_PRJNA2722046
#SBATCH --output=fastq_dump_PRJNA2722046_out

module load gcc/9.3.0
module load sra-toolkit/2.10.8

# for paired end reads only
fastq-dump --split-3  $1

# end
