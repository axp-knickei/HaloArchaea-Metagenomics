# Halo-Archaea Metagenomics Pipeline

**A High-Throughput Snakemake v8 Workflow for Hypersaline Ecosystems**

## 1. Overview and Biological Context
This pipeline is specialized for the metagenomic analysis of **Haloarchaea** (Class *Halobacteria*), which thrive in saturating salt concentrations (>2.5 M NaCl). Unlike standard bacterial pipelines, this workflow addresses unique halophilic adaptations:

* **"Salt-In" Strategy Verification:** Haloarchaea accumulate molar concentrations of KCl, resulting in a highly acidic proteome (median pI ~4.0–5.0). This pipeline includes a custom module to calculate proteome-wide isoelectric points (pI) to distinguish genuine Haloarchaea from "salt-out" bacteria or contaminants.
* **Gas Vesicle Preservation:** Standard quality control often removes biologically significant repetitive regions, such as gas vesicle gene clusters (*gvp*). This pipeline uses tuned `fastp` parameters to prevent the loss of these critical features.
* **Microdiversity Resolution:** Utilizes **MetaSPAdes** with multiple k-mer sizes to resolve strain heterogeneity common in halophilic communities.

## 2. Pipeline Architecture
The workflow is built on **Snakemake v8** using the **SLURM Executor Plugin** and a "Container-First" strategy for reproducibility.

| Stage | Tool | Notes |
| :--- | :--- | :--- |
| **QC** | `fastp` | Complexity threshold relaxed to 20 to preserve *gvp* genes. |
| **Assembly** | `MetaSPAdes` | Requires BigMem nodes (500GB+ RAM). |
| **Binning** | `SemiBin2` | Deep learning-based binning; accelerated via GPU. |
| **Taxonomy** | `GTDB-Tk` | v2.5.0+ using `skani` for ANI screening; optimized for Archaea. |
| **Quality** | `CheckM2` | Machine learning assessment (lineage-agnostic). |
| **Annotation**| `DRAM` | Metabolic profiling (KEGG/UniRef). |

## 3. Installation and Requirements

### Hardware Requirements
* **Assembly:** BigMem nodes (512 GB - 1 TB RAM).
* **Binning:** NVIDIA GPU (A100/V100 recommended).
* **Storage:** High-performance parallel file system (Lustre/GPFS) recommended for scratch space.

### Software Setup
The project uses `pyproject.toml` to manage core Python dependencies and individual `workflow/envs/*.yaml` files for tool-specific Conda environments.

1.  **Create the Master Environment** (Mamba recommended):
    ```bash
    mamba create -n snakemake_8 -c conda-forge -c bioconda \
        snakemake=8.4.6 \
        snakemake-executor-plugin-slurm \
        python=3.10 \
        pandas
    conda activate snakemake_8
    ```
    This environment provides the Snakemake orchestrator and general Python tools.

2.  **Install Project Core Dependencies**:
    The main CLI and helper scripts use dependencies defined in `pyproject.toml`.
    ```bash
    pip install .
    ```

3.  **Conda Environment Version Pinning**:
    All tool-specific Conda environments within `workflow/envs/` have their dependencies pinned to specific versions to ensure reproducibility. These environments will be automatically created/managed by Snakemake.
    **

2.  **Database Preparation**:
    * **GTDB-Tk:** Download R220 split package and extract.
    * **CheckM2:** Run `checkm2 database --download`.
    * **DRAM:** Run `DRAM-setup.py` (requires >100GB RAM).
    **

## 4. Configuration

### 1. Edit `config/config.yaml`
Define your project name, absolute paths to your downloaded databases, and key tool parameters.
This file now controls:
*   **Assembler selection:** Choose between `metaspades` and `megahit`.
*   **SemiBin2 environment:** Specify the binning model environment (e.g., `global`, `human_gut`).
*   **GTDB-Tk settings:** Configure ANI screening options (`use_skani`).
```yaml
databases:
  gtdb: "/path/to/gtdb_r220"
  checkm2: "/path/to/checkm2_db/uniref100.KO.1.dmnd"
  dram: "/path/to/DRAM_data"
assembler: "metaspades" # Options: metaspades, megahit
semibin2:
  environment: "global" # e.g., global, human_gut
gtdbtk:
  use_skani: True
```

### 2\. Edit `config/samples.tsv`

Map your sample names to their raw FASTQ files.

```tsv
sample	path_r1	path_r2
SaltPond_A	raw_data/A_R1.fastq.gz	raw_data/A_R2.fastq.gz
```

### 3\. Edit `config/resources.yaml`
This file is now loaded by the Snakemake workflow to manage computational resources (CPU, RAM, runtime, SLURM partitions). Match the partition names to your HPC's SLURM configuration.
```yaml
metaspades:
  slurm_partition: "bigmem" # Change to your cluster's high-memory partition
  mem_mb: 500000
```

## 5. Usage

Execute the pipeline using the provided `main.py` command-line interface (CLI). Ensure you have activated your master Snakemake environment and installed the project's core dependencies (`pip install .`). Apptainer/Singularity must be installed if you intend to use containers.

**Example: Local Execution**
```bash
python main.py --mode local --cores 8
```
This will run the workflow locally using 8 cores.

**Example: SLURM Execution**
```bash
python main.py --mode slurm --jobs 50 --config config/config.yaml
```
This will submit jobs to a SLURM cluster, allowing up to 50 concurrent jobs.

**Available `main.py` CLI options:**
```bash
python main.py --help
```

**Notes on SLURM and GPU:**
*   The `--singularity-args "--nv"` is automatically added in SLURM mode by the `main.py` CLI. This is critical for passing GPU access to tools like **SemiBin2** if they are run via a container.
*   The `config/resources.yaml` file defines default and rule-specific resources for SLURM job submission.

## 6. Output Structure

  * `results/assembly/`: Contigs and assembly graphs.
  * `results/binning/final_bins/`: Refined Metagenome-Assembled Genomes (MAGs).
  * `results/halo_analysis/`: **Proteomic pI Reports**. Check `*_proteome_pI.tsv` to verify the "salt-in" signature (median pI \< 5.0).
  * `results/annotation/dram/`: Metabolic heatmaps and gene inventories.
