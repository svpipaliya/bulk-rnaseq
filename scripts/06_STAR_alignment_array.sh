#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-user=shweta.pipaliya@epfl.ch
#SBATCH --array=1-101
#SBATCH --job-name=STAR_aln
#SBATCH --output=star_aln_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=star_aln_%A_%a.err # File to which STDERR will be written

# load all modules for runnning the STAR v.2.7 aligner
module load StdEnv/2020
module load star/2.7.9a

# print start date and echo hostname
start=`date +%s`
echo $HOSTNAME

# Define variables
GNM_DIR=/project/6001383/pip17/sars_cov_2/data/hg38_ref

# Path to the trimmed reads
FQ_DIR=/project/6001383/pip17/sars_cov_2/data/tools/salmon-1.9.0_linux_x86_64/bin

# Directory for output files
#BAM_DIR=/project/6001383/pip17/sars_cov_2/data/PRJEB43380/03_star_alignment

# set up job array parameters
sampleinfo="sample.sheet.txt"

# set variables to read the columns corresponding to sample names and r1 and r2
SAMPLE=`sed -n "$SLURM_ARRAY_TASK_ID"p $sampleinfo |  awk '{print $1}'`
r1="$SAMPLE"_trim_1.fq.gz
r2="$SAMPLE"_trim_2.fq.gz

# make output directory for results
#mkdir salmon_output_dir
#outdir= "salmon_output_dir"

#echo $SAMPLE

#if [ ! -e $outdir ]; then
    #mkdir $outdir
#fi

# parameter specifications for taking taking trimmed reads and performing alignment + quantification using STAR
STAR --genomeDir $GNM_DIR \
--runThreadN 1 \
--readFilesIn $r1 $r2 \
--readFilesCommand zcat \
--sjdbGTFfile $GNM_DIR/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFileNamePrefix $SAMPLE"_"

# print end date and echo total runtime
end=`date +%s`
runtime=$((end-start))
echo Runtime: $runtime seconds

# end