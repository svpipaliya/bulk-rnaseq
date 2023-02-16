
		
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
