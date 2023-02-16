rule rawFastqc:
	input:
		rawread="../input/rawReads/{fibro}_{frr}.fastq.gz"
	output:
		zip="../rawQC/{fibro}_{frr}_fastqc.gz", # the suffix _fastqc.zip is necessary for multiqc to find the file.
		html="../rawQC/{fibro}_{frr}_fastqc.html"
	params:
		path="../rawQC/"
	resources: 
		threads=1,
		runtime=240, 
		mem_mb=1024
	shell:
		"""
		fastqc {input.rawread} --threads {resources.threads} -o {params.path}
		"""
