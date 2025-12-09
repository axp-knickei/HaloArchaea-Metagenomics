# HaloArchaea Metagenomics Pipeline
**A High-Throughput Snakemake v8 Workflow for Hypersaline Ecosystems**

![Snakemake](https://img.shields.io/badge/snakemake-≥8.4.6-brightgreen.svg)
![Python](https://img.shields.io/badge/python-3.10-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green)
![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.17863847.svg)


![HaloArchaeLogo](HaloArchaea.png)

## Table of Contents
- [1. Overview](#1-overview)
- [2. Pipeline Architecture](#2-pipeline-architecture)  
- [3. Quick Start](#3-quick-start)
- [4. System Requirements](#4-system-requirements)
- [5. Installation & Setup](#5-installation--setup)
- [6. Configuration](#6-configuration)
- [7. Usage](#7-usage)
- [8. Output Structure](#8-output-structure)
- [9. Example Results](#9-example-results)
- [10. Troubleshooting](#10-troubleshooting)
- [11. References](#11-references)

## 1. Overview
This pipeline is specialized for the **metagenomic analysis of Haloarchaea** (Class *Halobacteria*), which thrive in saturating salt concentrations (>2.5 M NaCl). Hypersaline environments present unique challenges for bioinformatics that standard pipelines often miss.

**Key Features:**
* **"Salt-In" Strategy Verification:** Haloarchaea accumulate molar concentrations of KCl, resulting in a highly acidic proteome (median pI ~4.0–5.0). This workflow calculates proteome-wide isoelectric points to validate halophilic adaptations.
* **Gas Vesicle Preservation:** Standard QC often removes repetitive regions like gas vesicle gene clusters (*gvp*). We use tuned `fastp` parameters to preserve these critical features.
* **Microdiversity Resolution:** Uses **MetaSPAdes** with multiple k-mer sizes to resolve strain heterogeneity common in halophilic communities.
* **Modern Binning:** Implements **SemiBin2** with GPU acceleration and **CheckM2** for accurate quality assessment of Archaea.

## 2. Pipeline Architecture
The workflow is built on **Snakemake v8** using the **SLURM Executor Plugin** and a "Container-First" strategy for reproducibility.

| Stage | Tool | Description |
| :--- | :--- | :--- |
| **QC** | [**fastp**](https://github.com/OpenGene/fastp) | Quality control with relaxed complexity thresholds to preserve *gvp* genes. |
| **Assembly** | [**MetaSPAdes**](https://github.com/ablab/spades) | De novo assembly optimized for metagenomes (requires High RAM). |
| **Binning** | [**SemiBin2**](https://github.com/BigDataBiology/SemiBin) | Deep learning-based binning (GPU recommended). |
| **Taxonomy** | [**GTDB-Tk**](https://github.com/Ecogenomics/GTDBTk) | Taxonomic classification using the Genome Taxonomy Database (v2.5+). |
| **Quality** | [**CheckM2**](https://github.com/chklovski/CheckM2) | Machine learning-based completeness/contamination assessment. |
| **Annotation**| [**DRAM**](https://github.com/WrightonLabCSU/DRAM) | Distilled and Refined Annotation of Metabolism. |

![Workflow Halo Archaea](Pipeline_Halo_Archaea.png)

## 3. Quick Start
For users who want to run a test on a local machine immediately.

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/axp-knickei/HaloArchaea-Metagenomics.git](https://github.com/axp-knickei/HaloArchaea-Metagenomics.git)
    cd HaloArchaea-Metagenomics
    ```
2.  **Install the environment:**
    ```bash
    mamba create -n snakemake_8 -c conda-forge -c bioconda snakemake=8.4.6 snakemake-executor-plugin-slurm python=3.10 pandas
    mamba activate snakemake_8
    pip install .
    ```
3.  **Run the test:**
    ```bash
    # Runs the pipeline locally using 8 cores
    python main.py --mode local --cores 8
    ```

## 4. System Requirements

### Hardware Recommendations
Metagenomic assembly is computationally intensive.
* **Minimum (Testing/Small Data):** 64 GB RAM, 16 CPUs.
* **Production (Complex Metagenomes):**
    * **RAM:** 512 GB - 1 TB (Required for MetaSPAdes on large datasets).
    * **Storage:** 1 TB+ high-performance storage (Lustre/GPFS recommended).
    * **GPU:** NVIDIA GPU (A100/V100/T4) highly recommended for SemiBin2.

### Software Prerequisites
* **OS:** Linux (CentOS 7+, Ubuntu 20.04+).
* **Container Runtime:** Apptainer (formerly Singularity) is **required** for SLURM execution.
* **Package Manager:** [Mamba](https://github.com/mamba-org/mamba) or Conda.

## 🚀 5. Installation & Setup

### Step 1: Clone the Repository
```bash
git clone [https://github.com/axp-knickei/HaloArchaea-Metagenomics.git](https://github.com/axp-knickei/HaloArchaea-Metagenomics.git)
cd HaloArchaea-Metagenomics
````

### Step 2: Create the Master Environment

This environment runs the Snakemake orchestrator, not the heavy tools (which are handled via containers/envs).

```bash
mamba create -n snakemake_8 -c conda-forge -c bioconda \
    snakemake=8.4.6 \
    snakemake-executor-plugin-slurm \
    python=3.10 \
    pandas
mamba activate snakemake_8
```

### Step 3: Install Core Dependencies

```bash
pip install .
```

### Step 4: Download Required Databases

This pipeline relies on several large databases. Ensure you have \~200GB of disk space ready.

1.  **GTDB-Tk (R220):**
    ```bash
    wget [https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_r220_data.tar.gz](https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_r220_data.tar.gz)
    tar -xvzf gtdbtk_r220_data.tar.gz -C /path/to/databases/gtdb_r220
    ```
2.  **CheckM2:**
    ```bash
    checkm2 database --download --path /path/to/databases/checkm2_db
    ```
3.  **DRAM:**
    ```bash
    DRAM-setup.py prepare_databases --output_dir /path/to/databases/DRAM_data
    ```

## 6. Configuration

### 1\. Edit `config/samples.tsv`

Map your sample names to their raw FASTQ files.

```tsv
sample       path_r1                  path_r2
SaltPond_A   raw_data/A_R1.fastq.gz   raw_data/A_R2.fastq.gz
SaltPond_B   raw_data/B_R1.fastq.gz   raw_data/B_R2.fastq.gz
```

### 2\. Edit `config/config.yaml`

Set the paths to the databases you downloaded in Step 4.

```yaml
databases:
  gtdb: "/path/to/databases/gtdb_r220"
  checkm2: "/path/to/databases/checkm2_db/uniref100.KO.1.dmnd"
  dram: "/path/to/databases/DRAM_data"

assembler: "metaspades"      # Recommended for Haloarchaea
semibin2:
  environment: "global"      # Use 'global' for environmental samples
gtdbtk:
  use_skani: True            # Faster ANI calculation
```

## 7. Usage

### Local Execution (Small Test)

Run on a single machine without a scheduler.

```bash
python main.py --mode local --cores 8
```

### SLURM Cluster Execution (Production)

Submit jobs to a SLURM cluster. This automatically handles GPU requests for SemiBin2.

```bash
python main.py --mode slurm --jobs 50 --config config/config.yaml
```

**Note:** The pipeline automatically passes `--singularity-args "--nv"` to SemiBin2 jobs to enable GPU acceleration.

## 8. Output Structure

After the pipeline finishes, check the `results/` directory:

  * **`results/assembly/`**: Final contigs (`scaffolds.fasta`) and assembly graphs.
  * **`results/binning/final_bins/`**: High-quality MAGs (FASTA files).
  * **`results/halo_analysis/`**: `*_proteome_pI.tsv` (Crucial for validation).
  * **`results/taxonomy/`**: GTDB-Tk classification summaries.
  * **`results/annotation/dram/`**: Metabolic heatmaps (`product.html`) and gene inventories.

## 9. Example Results

### Proteome Isoelectric Point Validation

The file `results/halo_analysis/{sample}_proteome_pI.tsv` allows you to plot the pI distribution. A bimodal distribution with a peak \< 5.0 confirms the "Salt-In" strategy.

| protein\_id | pI | gravy | molecular\_weight |
| :--- | :--- | :--- | :--- |
| k141\_312\_1 | 4.23 | -0.15 | 24500 |
| k141\_312\_2 | 4.05 | -0.32 | 18200 |

### Metabolic Summary

Open `results/annotation/dram/{sample}/product.html` in your browser to view the heatmap of metabolic potential, including sulfur cycling and arsenic resistance genes common in hypersaline brines.

## 10. Troubleshooting

### Out of Memory (OOM) Errors

**Symptom:** Jobs fail with exit code 137 or "Out of memory" messages in logs.
**Solution:**

1.  **MetaSPAdes:** Increase `mem_gb` in `config/config.yaml` or `mem_mb` in `config/resources.yaml`.
2.  **Assembly:** If assembly still fails, try reducing the complexity by removing the largest k-mer (e.g., remove `127`) from the `kmer` list in `config/config.yaml`.

### GPU Not Detected

**Symptom:** SemiBin2 runs significantly slower or logs warnings about missing CUDA.
**Solution:**

1.  Ensure the `--nv` flag is being passed (handled automatically by `main.py` in SLURM mode).
2.  Verify your cluster supports Apptainer/Singularity with Nvidia support.
3.  Check `workflow/envs/binning.yaml` notes regarding PyTorch versions.

### Empty Bins

**Symptom:** The pipeline finishes, but `results/binning/final_bins/` is empty.
**Solution:**

  * Check `results/binning/semibin_output/` to see if SemiBin2 produced output.
  * Low biomass samples or poor assembly quality can result in zero high-quality bins. Check `results/qc/multiqc` to verify read quality.

## 11. References

1.  **Snakemake:** Mölder, F., et al. (2021). Sustainable data analysis with Snakemake. *F1000Research*, 10, 33.
2.  **fastp:** Chen, S., et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884–i890.
3.  **MetaSPAdes:** Nurk, S., et al. (2017). metaSPAdes: a new versatile metagenomic assembler. *Genome Research*, 27(5), 824–834.
4.  **SemiBin2:** Pan, S., et al. (2023). SemiBin2: self-supervised contrastive learning leads to better MAGs. *Bioinformatics*, 39(Suppl 1), i21–i29.
5.  **GTDB-Tk:** Chaumeil, P. A., et al. (2022). GTDB-Tk v2: memory friendly classification with the Genome Taxonomy Database. *Bioinformatics*, 38(23), 5315–5316.
6.  **CheckM2:** Chklovski, A., et al. (2023). CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nature Methods*, 20, 1203–1212.
7.  **DRAM:** Shaffer, M., et al. (2020). DRAM for distilling microbial metabolism to automate the curation of microbiome function. *Nucleic Acids Research*, 48(16), 8883–8900.

## 12. Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{HaloArcaea-Metagenomics,
  author = {Alex Prima},
  title = {HaloArchaea-Metagenomics Pipeline, a high-throughput Snakemake v8 workflow designed for hypersaline ecosystems},
  year = {2025},
  url = {https://github.com/axp-knickei/HaloArchaea-Metagenomics},
  doi = {10.5281/zenodo.17863847}
}
```
Also cite the individual tools (see References section above).