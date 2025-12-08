rule metaspades:
    input:
        r1 = "results/qc/fastp/{sample}_R1.trimmed.fastq.gz",
        r2 = "results/qc/fastp/{sample}_R2.trimmed.fastq.gz"
    output:
        contigs = "results/assembly/{sample}/contigs.fasta",
        graph = "results/assembly/{sample}/assembly_graph.fastg"
    log:
        "logs/metaspades/{sample}.log"
    threads: 24
    resources:
        mem_mb = 500000, 
        runtime = 2880
    conda:
        "../envs/assembly.yaml"
    params:
        k = ",".join(map(str, config["metaspades"]["kmer"])),
        outdir = "results/assembly/{sample}"
    shell:
        """
        metaspades.py \
            -1 {input.r1} -2 {input.r2} \
            -o {params.outdir} \
            -k {params.k} \
            -t {threads} \
            --memory 490 \
            > {log} 2>&1
        
        mv {params.outdir}/contigs.fasta {output.contigs}
        mv {params.outdir}/assembly_graph.fastg {output.graph}
        """

rule index_contigs:
    input:
        "results/assembly/{sample}/contigs.fasta"
    output:
        "results/assembly/{sample}/contigs.fasta.1.bt2"
    threads: 8
    conda:
        "../envs/binning.yaml"
    shell:
        "bowtie2-build {input} results/assembly/{sample}/contigs.fasta --threads {threads}"