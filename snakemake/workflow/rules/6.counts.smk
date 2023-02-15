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
			
# this rule is only applicable if you want to perform transcript-level quantification using Salmon
#rule salmonQuant:
	#input:
		#bam=expand("../output/salmonAligned/{fibro}.Salmon.Aligned.bam", fibro = fibro),
		#ref="../resources/hg38/hg38_transcriptome.fa"
	#output:
		#salmonQuant=expand("../output/salmonQuant/{fibro}_quant.sf", fibro=fibro)
	#resources:
		#threads=5,
		#runtime=4320, 
		#mem_mb=4096
	#shell:
		#"""
		#salmon quant -t {input.ref} -l A \
		#-a {input.bam} -o {output.salmonQuant}
		#"""
