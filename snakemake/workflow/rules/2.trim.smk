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
		
rule fastp:
	input:
		read1="../input/rawReads/{fibro}_1.fastq.gz",
		read2="../input/rawReads/{fibro}_2.fastq.gz"
	output:
		read1T="../output/trimmedReads/{fibro}_trim_1.fastq.gz",
		read2T="../output/trimmedReads/{fibro}_trim_2.fastq.gz",
		html="../output/trimmedReads/{fibro}_fastq.html"
	resources: 
		threads=4,
		runtime=2880, 
		mem_mb=4096
	shell:
		"""
		fastp -i {input.read1} -I {input.read2} \
		-o {output.read1T} -O {output.read2T} \
		 -h {output.html} --detect_adapter_for_pe -g -x -p 
		"""
