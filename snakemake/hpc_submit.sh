#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=3850
#SBATCH --time=72:00:00
#SBATCH --output=log/main_%j

###### activate conda environment ######
eval "$(conda shell.bash hook)"
conda activate rnaseq

##### Run the pipaline #####
snakemake --profile ./env/slurm --keep-going --rerun-incomplete --nolock --latency-wait 999

##### end #####