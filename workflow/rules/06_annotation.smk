rule dram_annotate:
    input:
        bins = "results/binning/final_bins/{sample}"
    output:
        annotations = "results/annotation/dram/{sample}/annotations.tsv",
        genes = "results/annotation/dram/{sample}/genes.faa"
    threads: 24
    resources:
        mem_mb = 100000
    container:
        "docker://wrightonlab/dram:latest"
    params:
        db_dir = config["databases"]["dram"],
        out_dir = "results/annotation/dram/{sample}"
    shell:
        """
        DRAM.py annotate \
            --input_fasta_dir {input.bins} \
            --output_dir {params.out_dir} \
            --config_loc {params.db_dir}/DRAM_config.tsv \
            --threads {threads}
        """