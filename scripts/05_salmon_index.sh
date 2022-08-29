#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH	--cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=shweta.pipaliya@epfl.ch
#SBATCH --job-name=Salmon_transcriptome_index
#SBATCH --output=Salmon__index_out
#SBATCH --error=Salmon_index_err

# set path to salmon executable script
set SM = /work/backup/gr-fe/pipaliya/tools/salmon-1.9.0_linux_x86_64/bin/salmon

# parameter specification for running STAR to generate a genome index
$SM index -t hg38_transcriptome.fa.gz -i hg38_index

#end