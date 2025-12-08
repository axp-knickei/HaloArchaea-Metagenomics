# Installation & System Requirements

## Hardware Recommendations
Metagenomic assembly of hypersaline communities is computationally intensive due to the complexity of the data.

* **Minimum (Testing/Small Data):** * 64 GB RAM
    * 16 CPUs
* **Production (Complex Metagenomes):**
    * **RAM:** 512 GB - 1 TB (Required for MetaSPAdes on large datasets).
    * **Storage:** 1 TB+ high-performance storage (Lustre/GPFS recommended).
    * **GPU:** NVIDIA GPU (A100/V100/T4) is **highly recommended** for SemiBin2.

## Software Prerequisites
* **OS:** Linux (CentOS 7+, Ubuntu 20.04+).
* **Container Runtime:** [Apptainer](https://apptainer.org/) (formerly Singularity) is **required** for SLURM execution.
* **Package Manager:** [Mamba](https://github.com/mamba-org/mamba) or Conda.

## Step-by-Step Installation

### 1. Clone the Repository
```bash
git clone [https://github.com/axp-knickei/HaloArchaea-Metagenomics.git](https://github.com/axp-knickei/HaloArchaea-Metagenomics.git)
cd HaloArchaea-Metagenomics
````

### 2\. Create the Master Environment

This environment runs the Snakemake orchestrator. Heavy tools (SPAdes, GTDB-Tk) run in separate containers/environments automatically.

```bash
mamba create -n snakemake_8 -c conda-forge -c bioconda \
    snakemake=8.4.6 \
    snakemake-executor-plugin-slurm \
    python=3.10 \
    pandas

mamba activate snakemake_8
```

### 3\. Install Core Dependencies

```bash
pip install .
```