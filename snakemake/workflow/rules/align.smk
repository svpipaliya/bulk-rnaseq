rule align:
    input:
        fq1="reads/{sample}_R1.1.fastq",
        fq2="reads/{sample}_R2.1.fastq",
        index="resources/star_genome",
        idx="resources/genome.gtf"
    output:
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        "results/star/{sample}-{unit}/ReadsPerGene.out.tab",
    log:
        "logs/star/{sample}-{unit}.log",
    params:
        index=lambda wc, input: input.index,
        extra="--readFilesCommand zcat --outSAMtype BAM Unsorted SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --quantMode GeneCounts --sjdbGTFfile {} {}".format(
            "resources/genome.gtf", config["params"]["star"]
        ),
    threads: 8
    wrapper:
        "v1.21.1/bio/star/align"