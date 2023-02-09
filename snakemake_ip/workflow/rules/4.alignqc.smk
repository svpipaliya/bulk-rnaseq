# generate reports of the final qc'd file
rule multiqc:
	input:
		bam=expand("../output/starAligned/{fibro}.Aligned.sortedbyCoord.out.bam", fibro=fibro)
	output:
		qc="../output/bamQC/multiqc_star_report.html"
	resources:
		threads=1,
		runtime=60,
		mem_mb=1024
	shell:
		"""
		multiqc --force {input.bam} -o {output.qc}
		"""	