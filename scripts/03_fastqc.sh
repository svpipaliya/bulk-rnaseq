#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH	--cpus-per-task=1
#SBATCH --mem-per-cpu=20000M
#SBATCH --mail-user=shweta.pipaliya@epfl.ch
#SBATCH --job-name=fastqc_GSE172114
#SBATCH --output=fastqc_GSE172114_out

# Loading modules required for script commands
module load nixpkgs/16.09
module load fastqc/0.11.9

#  Running FASTQC
fastqc -t 6 *.fastq.gz


# make a results directory
mkdir /results/fastqc

# Moving files to our results directory
mv *fastqc* ../results/fastqc/

# end
