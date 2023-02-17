rule multiqc:
	input:
		bam=expand("../starAligned/{fibro}.Aligned.sortedbyCoord.out.bam", fibro=fibro)
	output:
		qc="../multiQC/multiqc_star_report.html"
	resources:
		threads=1,
		runtime=60,
		mem_mb=1024
	shell:
		"""
		multiqc --force -o "../multiQC/" -n {output.qc}
		"""	
