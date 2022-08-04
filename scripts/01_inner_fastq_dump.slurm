#!/bin/bash
#SBATCH --time=5-00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH	--cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-user=shweta.pipaliya@epfl.ch
#SBATCH --job-name=fastq_dump_GSE172114
#SBATCH --output=fastq_dump_GSE172114_out

module load nixpkgs/16.09
module load sra-toolkit/2.9.6

# initialize a variable with an intuitive name to store the name of the input fastq file

fq=$1
fastq $fq

# grab base of filename for naming outputs

base=`basename $fq .subset.fq`
echo "Sample name is $base" 

# for paired end reads only
fastq-dump --split-3  $1

# zip all resulting read files
gzip *.fastq

# end


