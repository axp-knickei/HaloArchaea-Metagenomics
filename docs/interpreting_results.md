# Interpreting Results

After a successful run, the `results/` directory will contain the following:

## 1. Quality Control (`results/qc/`)
* **`multiqc/multiqc_report.html`**: Aggregate report of read quality, adapter trimming, and duplication rates. Check this first to ensure raw data quality was sufficient.

## 2. Assembly (`results/assembly/`)
* **`{sample}/contigs.fasta`**: The final assembled contigs.
* **`{sample}/assembly_graph.fastg`**: The assembly graph, useful for visualizing strain heterogeneity in tools like Bandage.

## 3. Bins (`results/binning/`)
* **`final_bins/{sample}/*.fa`**: High-quality Metagenome-Assembled Genomes (MAGs). Each fasta file represents a distinct genome bin.

## 4. Taxonomy & Quality (`results/taxonomy/`, `results/quality/`)
* **`taxonomy/gtdbtk/{sample}/classify/gtdbtk.bac120.summary.tsv`**: The taxonomic identity of each bin. Look for Class *Halobacteria*.
* **`quality/checkm2/{sample}/quality_report.tsv`**: Look for:
    * `Completeness`: > 90% is excellent.
    * `Contamination`: < 5% is excellent.

## 5. Haloarchaeal Analysis (`results/halo_analysis/`)
* **`{sample}_proteome_pI.tsv`**: A tab-separated file listing the pI of every protein.
    * **Usage:** Plot a histogram of the `pI` column. A bimodal distribution with a major peak around 4.5 indicates a salt-in adapted community.

## 6. Functional Annotation (`results/annotation/`)
* **`dram/{sample}/product.html`**: An interactive heatmap of metabolic functions.
    * **Look for:** Sulfur cycling genes, Arsenic resistance, and Rhodopsins (light-driven proton pumps) common in haloarchaea.