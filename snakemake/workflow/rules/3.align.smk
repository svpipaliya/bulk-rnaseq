rule starAlign:
	input:
		read1="../output/trimmedReads/{fibro}_trim_1.fastq.gz",
		read2="../output/trimmedReads/{fibro}_trim_2.fastq.gz",
		ref="hg38/",
		annotation="hg38/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf"
	output:
		bam="../output/starAligned/{fibro}.Aligned.sortedbyCoord.out.bam",
		log="../output/starAligned/{fibro}.Log.final.out"
	params: 
		prefix="starAligned/{fibro}"
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
		--outReadsUnmapped Within --outSAMattributes Standard 
		"""

# this rule is only applicable if you want to generate alignments using Salmon
#rule salmonAlign:
	#input:
		#read1=expand("../output/trimmedReads/{fibro}_trim_1.fastq.gz",fibro=fibro),
		#read2=expand("../output/trimmedReads/{fibro}_trim_2.fastq.gz",fibro=fibro),
		#index="salmonIndex/salmon_hg38_index"
	#output:
		#salmonAligned=expand("../output/salmonAligned/{fibro}.Salmon.Aligned.bam", fibro = fibro)
	#resources: 
		#threads=5,
		#runtime=4320, 
		#mem_mb=8192
	#shell:
		#"""
		#salmon quant -i {input.index} -l A -p 8 \
		#-1 {input.read1} -2 {input.read2} --validateMappings \
		#--seqBias --useVBOpt -o {output.salmonAligned} 
		#"""
