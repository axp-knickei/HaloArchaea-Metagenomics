rule extract_genes:
    input:
        "results/annotation/dram/{sample}/genes.faa"
    output:
        "results/halo_analysis/{sample}_genes.faa"
    shell:
        "cp {input} {output}"

rule calculate_proteome_pi:
    input:
        "results/halo_analysis/{sample}_genes.faa"
    output:
        "results/halo_analysis/{sample}_proteome_pI.tsv"
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/calculate_pI.py"