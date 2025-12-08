rule gtdbtk_classify:
    input:
        bins_dir = "results/binning/final_bins/{sample}"
    output:
        summary = "results/taxonomy/gtdbtk/{sample}/classify/gtdbtk.bac120.summary.tsv"
    threads: 24
    resources:
        mem_mb = 120000
    container:
        "docker://ecogenomic/gtdbtk:2.5.0"
    params:
        outdir = "results/taxonomy/gtdbtk/{sample}",
        db = config["databases"]["gtdb"]
    shell:
        """
        export GTDBTK_DATA_PATH={params.db}
        gtdbtk classify_wf \
            --genome_dir {input.bins_dir} \
            --out_dir {params.outdir} \
            --extension fa \
            --cpus {threads} \
            --skip_ani_screen 
        """