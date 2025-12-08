# Assembly Rules
# Supports MetaSPAdes and MEGAHIT based on config["assembler"]

if config["assembler"] == "metaspades":
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
            mem_mb = config["metaspades"]["mem_gb"] * 1000,
            runtime = 2880
        conda:
            "../envs/assembly.yaml"
        params:
            k = ",".join(map(str, config["metaspades"]["kmer"])),
            outdir = "results/assembly/{sample}",
            mem = config["metaspades"]["mem_gb"]
        shell:
            """
            metaspades.py \
                -1 {input.r1} -2 {input.r2} \
                -o {params.outdir} \
                -k {params.k} \
                -t {threads} \
                --memory {params.mem} \
                > {log} 2>&1
            
            # Ensure output matches expectation
            if [ -f {params.outdir}/scaffolds.fasta ];
            then
                mv {params.outdir}/scaffolds.fasta {output.contigs}
            else
                echo "Error: Scaffolds not found" >&2
                exit 1
            fi
            
            # Graph might not exist if assembly failed early but contigs exist (unlikely but possible)
            if [ -f {params.outdir}/assembly_graph.fastg ]; 
            then
                touch {output.graph}
            fi
            """

elif config["assembler"] == "megahit":
    rule megahit:
        input:
            r1 = "results/qc/fastp/{sample}_R1.trimmed.fastq.gz",
            r2 = "results/qc/fastp/{sample}_R2.trimmed.fastq.gz"
        output:
            contigs = "results/assembly/{sample}/contigs.fasta",
            graph = "results/assembly/{sample}/assembly_graph.fastg" # Mocking graph for compatibility
        log:
            "logs/megahit/{sample}.log"
        threads: 24
        resources:
            mem_mb = 64000, # Default reasonable value
            runtime = 1440
        conda:
            "../envs/assembly.yaml" # Assuming megahit is in the same env or I should create a new one.
                                   # ideally separate, but for now assuming user manages envs.
        params:
            outdir = "results/assembly/{sample}",
            out_prefix = "megahit"
        shell:
            """
            rm -rf {params.outdir} # Megahit requires empty dir
            megahit \
                -1 {input.r1} -2 {input.r2} \
                -o {params.outdir} \
                -t {threads} \
                --out-prefix {params.out_prefix} \
                > {log} 2>&1
            
            mv {params.outdir}/{params.out_prefix}.contigs.fa {output.contigs}
            touch {output.graph} # Megahit doesn't produce FASTG by default in the same way
            """

else:
    raise ValueError(f"Unknown assembler: {config['assembler']}")

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
