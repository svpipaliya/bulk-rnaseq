rule fastp_pe:
    input:
        sample=["reads/pe/{sample}.1.fastq", "reads/pe/{sample}.2.fastq"]
    output:
        trimmed=["trimmed/pe/{sample}.1.fastq", "trimmed/pe/{sample}.2.fastq"],
        html="report/pe/{sample}.html"
    log:
        "logs/fastp/pe/{sample}.log"
    params:
        adapters="--detect_adapter_for_pe",
        extra="--trim_poly_g --trim_poly_x --overrepresentation_analysis"
    threads: 2
    wrapper:
        "v1.21.1/bio/fastp"