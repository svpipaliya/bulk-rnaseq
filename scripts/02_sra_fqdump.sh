#!/bin/bash

# for every SRR in the list of SRRs file
for srr in $(cat SRR_Acc_List.txt)
  do
  # call the bash script that does the fastq dump, passing it the SRR number next in file
  sbatch 01_inner_fastq_dump.slurm $srr
  sleep 1	# wait 1 second between each job submission
done