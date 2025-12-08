rule map_reads:
    input:
        idx = "results/assembly/{sample}/contigs.fasta.1.bt2",
        contigs = "results/assembly/{sample}/contigs.fasta",
        r1 = "results/qc/fastp/{sample}_R1.trimmed.fastq.gz",
        r2 = "results/qc/fastp/{sample}_R2.trimmed.fastq.gz"
    output:
        bam = temp("results/mapping/{sample}.sorted.bam")
    threads: 16
    conda:
        "../envs/binning.yaml"
    shell:
        """
        bowtie2 -x {input.contigs} \
            -1 {input.r1} -2 {input.r2} \
            -p {threads} | \
        samtools view -bS - | \
        samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """

rule semibin2_binning:
    input:
        contigs = "results/assembly/{sample}/contigs.fasta",
        bam = "results/mapping/{sample}.sorted.bam"
    output:
        bins_dir = directory("results/binning/final_bins/{sample}")
    threads: 16
    resources:
        gpu = 1,
        mem_mb = 32000
    conda:
        "../envs/binning.yaml"
    params:
        environment = config["semibin2"]["environment"],
        tmp_out = "results/binning/semibin_output/{sample}"
    shell:
        """
        SemiBin2 single_easy_bin \
            -i {input.contigs} \
            -b {input.bam} \
            -o {params.tmp_out} \
            --environment {params.environment} \
            --threads {threads}
        
        mkdir -p {output.bins_dir}

        # Check if the output directory exists and is not empty
        if [ -d "{params.tmp_out}/output_bins" ]; then
            # Copy only if files exist, suppressing errors if no match
            cp {params.tmp_out}/output_bins/*.fa {output.bins_dir}/ 2>/dev/null || true
            
            # Optional: Check if we actually copied anything
            if [ -z "$(ls -A {output.bins_dir})" ]; then
                 echo "Warning: SemiBin2 finished but produced no bins for {wildcards.sample}" >&2
            fi
        else
            echo "Warning: SemiBin2 output directory not found for {wildcards.sample}" >&2
        fi