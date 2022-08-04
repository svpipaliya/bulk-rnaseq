#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH	--cpus-per-task=16
#SBATCH --mem-per-cpu=2000M
#SBATCH --mail-user=shweta.pipaliya@epfl.ch
#SBATCH --job-name=md5sums_pe_reads
#SBATCH --output=md5sums_pe_reads_out

#  Run md5sums on every output and parse out hashes in a text file
md5sum *.gz

# end