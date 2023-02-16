rule featurecounts:
	input:
		bam=expand("../output/starAligned/{fibro}.Aligned.sortedbyName.out.bam", fibro=fibro),
		annotation="../resources/hg38/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf"
	output:
		"../output/featureCounts/hgid_feature_counts.txt",
	resources:
		threads=8, 
		runtime=4320, 
		mem_mb=4096
	shell:
		"""
		featureCounts {input.bam} -a {input.annotation} -s 2 -p \
		--countReadPairs -t exon -T {resources.threads} -g gene_id -o {output}
		"""
			
