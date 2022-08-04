#!/bin/bash
#SBATCH --time=5-00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-user=shweta.pipaliya@epfl.ch
#SBATCH --job-name=fastp_GSE172114
#SBATCH --output=fastp_GSE172114_out

# Loading modules required for script commands
module load StdEnv/2020
module load fastp/0.23.1

# for every SRR in the list of SRRs file
for srr in $(cat SRR_Acc_List.txt)
  do
  echo -e $srr
 
     fastp -i ${srr}_1.fastq.gz\
     -I ${srr}_2.fastq.gz\
     -o ${srr}_trim_1.fastq.gz\
     -O ${srr}_trim_2.fastq.gz\
     -h ${srr}_report.html.gz\
     -j ${srr}_report.json.gz\
     --detect_adapter_for_pe\
     --thread 1\
     -g -x -p
done