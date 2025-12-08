rule checkm2_quality:
    input:
        bins = "results/binning/final_bins/{sample}"
    output:
        summary = "results/quality/checkm2/{sample}/quality_report.tsv"
    log:
        "logs/checkm2/{sample}.log"
    threads: 24
    resources:
        mem_mb = 120000
    container:
        "docker://chklovski/checkm2:latest"
    params:
        db = config["databases"]["checkm2"],
        outdir = "results/quality/checkm2/{sample}"
    shell:
        """
        checkm2 predict \
            --threads {threads} \
            --input {input.bins} \
            --output-directory {params.outdir} \
            --database_path {params.db} \
            --force
        """