#!/bin/bash
#SBATCH --time=7-00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH	--cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-user=shweta.pipaliya@epfl.ch
#SBATCH --job-name=fastqc_GSE172114
#SBATCH --output=fastqc_GSE172114_out

# Loading modules required for script commands
module load nixpkgs/16.09
module load fastqc/0.11.9

#  Running FASTQC
fastqc -t 6 *.fastq.gz

# end
