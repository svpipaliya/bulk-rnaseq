#!/bin/bash
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem-per-cpu=20G
#SBATCH --mail-user=shweta.pipaliya@epfl.ch
#SBATCH --array=1-101
#SBATCH --job-name=Salmon_aln_quant
#SBATCH --output=slurmout/salmon_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/salmon_%A_%a.err # File to which STDERR will be written

# print start date and echo hostname
start=`date +%s`
echo $HOSTNAME

# set up job array parameters
sampleinfo="sample.sheet.txt"

# set variables to read the columns corresponding to sample names and r1 and r2
SAMPLE=`sed -n "$SLURM_ARRAY_TASK_ID"p $sampleinfo |  awk '{print $1}'`
r1=`sed -n "$SLURM_ARRAY_TASK_ID"p $sampleinfo |  awk '{print $2}'`
r2=`sed -n "$SLURM_ARRAY_TASK_ID"p $sampleinfo |  awk '{print $3}'`

# make output directory for results
mkdir salmon_output_dir
outdir= "salmon_output_dir"

echo $SAMPLE

if [ ! -e $outdir ]; then
    mkdir $outdir
fi
 
# parameter specifications for taking taking trimmed reads and performing alignment + quantification using Salmon
./salmon quant -i hg38_index \
-l A \
-p 8 \
-1 $r1 \
-2 $r2 \
--validateMappings \
--seqBias \
--useVBOpt \
-o $outdir/$SAMPLE.out

# print end date and echo total runtime
end=`date +%s`
runtime=$((end-start))
echo Runtime: $runtime seconds

# end