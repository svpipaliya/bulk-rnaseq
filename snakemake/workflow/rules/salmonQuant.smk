rule salmonQuant:
	input:
		bam=expand("../output/salmonAligned/{fibro}.Salmon.Aligned.bam", fibro = fibro),
		ref="../resources/hg38/hg38_transcriptome.fa"
	output:
		salmonQuant=expand("../output/salmonQuant/{fibro}_quant.sf", fibro=fibro)
	resources:
		threads=5,
		runtime=4320, 
		mem_mb=4096
	shell:
		"""
		salmon quant -t {input.ref} -l A \
		-a {input.bam} -o {output.salmonQuant}
		"""
