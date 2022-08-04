#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20000M
#SBATCH --mail-user=shweta.pipaliya@epfl.ch
#SBATCH --job-name=cp_pe_reads_GSE172114
#SBATCH --output=cp_pe_reads_GSE172114

cp -a /home/pip17/projects/def-dacks/pip17/sars_cov_2/data/GSE172114/pe_reads/. /home/pip17/scratch/sars_cov_2/data/GSE172114/fastqc_qc