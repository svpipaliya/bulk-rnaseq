rule salmonAlign:
	input:
		read1=expand("../trimmedReads/{fibro}_trim_1.fastq.gz",fibro=fibro),
		read2=expand("../trimmedReads/{fibro}_trim_2.fastq.gz",fibro=fibro),
		index="salmonIndex/salmon_hg38_index"
	output:
		salmonAligned=expand("../salmonAligned/{fibro}.Salmon.Aligned.bam", fibro = fibro)
	resources: 
		threads=5,
		runtime=4320, 
		mem_mb=8192
	shell:
		"""
		salmon quant -i {input.index} -l A -p 8 \
		-1 {input.read1} -2 {input.read2} --validateMappings \
		--seqBias --useVBOpt -o {output.salmonAligned} 
		"""
