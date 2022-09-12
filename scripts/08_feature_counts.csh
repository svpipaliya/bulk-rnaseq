#!/bin/csh
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH	--cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-user=shweta.pipaliya@epfl.ch
#SBATCH --job-name=featureCounts_GSE172114
#SBATCH --output=featureCounts_GSE172114.o
#SBATCH --error=featureCounts_GSE172114.e

## load modules
module load subread/2.0.3
module load StdEnv/2020 
module load gcc/9.3.0

## set other variables 
set REF_GTF = /project/6001383/pip17/sars_cov_2/data/GSE172114/hg38_ref/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf
#set featureCounts = /project/6001383/pip17/sars_cov_2/data/GSE172114/tools/subread-2.0.3-Linux-x86_64/bin/featureCounts

## make output directory
#mkdir fc_GSE172114.output

## run FeatureCounts on the BAM files to generate gene counts for reads
featureCounts -a $REF_GTF \
 -s 2 \
 -p --countReadPairs \
 -t exon \
 -T 8 \
 -g gene_id \
 -o GSE172114_counts.txt \
 SRR14240730_mapAligned.sortedByCoord.out.bam SRR14240731_mapAligned.sortedByCoord.out.bam SRR14240732_mapAligned.sortedByCoord.out.bam \
 SRR14240733_mapAligned.sortedByCoord.out.bam SRR14240734_mapAligned.sortedByCoord.out.bam SRR14240735_mapAligned.sortedByCoord.out.bam \
 SRR14240736_mapAligned.sortedByCoord.out.bam SRR14240737_mapAligned.sortedByCoord.out.bam SRR14240738_mapAligned.sortedByCoord.out.bam \
 SRR14240740_mapAligned.sortedByCoord.out.bam SRR14240741_mapAligned.sortedByCoord.out.bam SRR14240742_mapAligned.sortedByCoord.out.bam \
 SRR14240744_mapAligned.sortedByCoord.out.bam SRR14240745_mapAligned.sortedByCoord.out.bam SRR14240746_mapAligned.sortedByCoord.out.bam \
 SRR14240747_mapAligned.sortedByCoord.out.bam SRR14240748_mapAligned.sortedByCoord.out.bam SRR14240749_mapAligned.sortedByCoord.out.bam \
 SRR14240750_mapAligned.sortedByCoord.out.bam SRR14240751_mapAligned.sortedByCoord.out.bam SRR14240752_mapAligned.sortedByCoord.out.bam \
 SRR14240753_mapAligned.sortedByCoord.out.bam SRR14240754_mapAligned.sortedByCoord.out.bam SRR14240755_mapAligned.sortedByCoord.out.bam \
 SRR14240756_mapAligned.sortedByCoord.out.bam SRR14240757_mapAligned.sortedByCoord.out.bam SRR14240758_mapAligned.sortedByCoord.out.bam \
 SRR14240759_mapAligned.sortedByCoord.out.bam SRR14240760_mapAligned.sortedByCoord.out.bam SRR14240761_mapAligned.sortedByCoord.out.bam \
 SRR14240762_mapAligned.sortedByCoord.out.bam SRR14240763_mapAligned.sortedByCoord.out.bam SRR14240764_mapAligned.sortedByCoord.out.bam \
 SRR14240765_mapAligned.sortedByCoord.out.bam SRR14240766_mapAligned.sortedByCoord.out.bam SRR14240767_mapAligned.sortedByCoord.out.bam \
 SRR14240768_mapAligned.sortedByCoord.out.bam SRR14240769_mapAligned.sortedByCoord.out.bam SRR14240770_mapAligned.sortedByCoord.out.bam \
 SRR14240771_mapAligned.sortedByCoord.out.bam SRR14240778_mapAligned.sortedByCoord.out.bam SRR14240774_mapAligned.sortedByCoord.out.bam \
 SRR14240779_mapAligned.sortedByCoord.out.bam SRR14240780_mapAligned.sortedByCoord.out.bam SRR14240782_mapAligned.sortedByCoord.out.bam \
 SRR14240783_mapAligned.sortedByCoord.out.bam SRR14240785_mapAligned.sortedByCoord.out.bam SRR14240775_mapAligned.sortedByCoord.out.bam \
 SRR14240786_mapAligned.sortedByCoord.out.bam SRR14240787_mapAligned.sortedByCoord.out.bam SRR14240788_mapAligned.sortedByCoord.out.bam \
 SRR14240789_mapAligned.sortedByCoord.out.bam SRR14240790_mapAligned.sortedByCoord.out.bam SRR14240791_mapAligned.sortedByCoord.out.bam \
 SRR14240792_mapAligned.sortedByCoord.out.bam SRR14240793_mapAligned.sortedByCoord.out.bam SRR14240794_mapAligned.sortedByCoord.out.bam \
 SRR14240795_mapAligned.sortedByCoord.out.bam SRR14240796_mapAligned.sortedByCoord.out.bam SRR14240797_mapAligned.sortedByCoord.out.bam \
 SRR14240798_mapAligned.sortedByCoord.out.bam SRR14240772_mapAligned.sortedByCoord.out.bam SRR14240773_mapAligned.sortedByCoord.out.bam \
 SRR14240776_mapAligned.sortedByCoord.out.bam SRR14240777_mapAligned.sortedByCoord.out.bam \
  
## end



