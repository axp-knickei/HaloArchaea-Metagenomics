rule fastp_pe:
    input:
        r1 = "raw_data/{sample}_R1.fastq.gz",
        r2 = "raw_data/{sample}_R2.fastq.gz"
    output:
        r1 = "results/qc/fastp/{sample}_R1.trimmed.fastq.gz",
        r2 = "results/qc/fastp/{sample}_R2.trimmed.fastq.gz",
        json = "results/qc/fastp/{sample}.json",
        html = "results/qc/fastp/{sample}.html"
    log:
        "logs/fastp/{sample}.log"
    threads: 8
    resources:
        mem_mb = 8000,
        runtime = 60
    conda:
        "../envs/qc.yaml"
    params:
        qualified_quality = config["fastp"]["qualified_quality_phred"],
        length_req = config["fastp"]["length_required"],
        complexity = config["fastp"]["complexity_threshold"],
        extra = config["fastp"]["extra"]
    shell:
        """
        fastp \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --qualified_quality_phred {params.qualified_quality} \
            --length_required {params.length_req} \
            --low_complexity_filter \
            --complexity_threshold {params.complexity} \
            {params.extra} \
            --json {output.json} --html {output.html} \
            --thread {threads} \
            2> {log}
        """

rule multiqc:
    input:
        expand("results/qc/fastp/{sample}.json", sample=SAMPLES)
    output:
        "results/qc/multiqc/multiqc_report.html"
    conda:
        "../envs/qc.yaml"
    shell:
        "multiqc results/qc/fastp -o results/qc/multiqc"