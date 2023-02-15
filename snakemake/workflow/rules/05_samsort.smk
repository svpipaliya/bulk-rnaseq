rule bamsort:
	input:
		bam=rules.starAlign.output.bam
	output:
		"../output/samsort/{fibro}.Aligned.sortedByName.out.bam"
	resources:
		threads=8,
		runtime=1440,
		mem_mb=1024
	shell:
		"""
		samtools sort -n @ {resources.threads} -f {input.bam} -o {output}
		"""
