rule salmonIndex:
	input:
		ref="../resources/hg38/hg38_transcriptome.fa"
	output:
		directory("../resources/salmon_hg38_index")
	resources: 
		threads=8,
		runtime=1440, 
		mem_mb=8192
	shell:
		"""
		salmon index -t {input.ref} -i {output}
		"""
