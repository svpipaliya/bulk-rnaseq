rule starAlign:
	input:
		read1="../output/trimmedReads/{fibro}_trim_1.fastq.gz",
		read2="../output/trimmedReads/{fibro}_trim_2.fastq.gz",
		ref="hg38/",
		annotation="hg38/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf"
	output:
		bam="../output/starAligned/{fibro}.Aligned.sortedbyCoord.out.bam",
		log="../output/starAligned/{fibro}.Log.final.out"
	log:
		"../log/starAligned/{fibro}.log"
	params: 
		prefix="../output/starAligned/{fibro}"
	resources: 
		threads=1, 
		runtime=4320, 
		mem_mb=64000
	shell:
		"""
		STAR --runThreadN 1 --genomeDir {input.ref} \
		--readFilesIn {input.read1} {input.read2} \
		--readFilesCommand zcat --sjdbGTFfile {input.annotation} \
		--outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate \
		--limitBAMsortRAM 16000000000 \
		--outReadsUnmapped Within --outSAMattributes Standard > {output.bam} {output.log} {log}
		"""

