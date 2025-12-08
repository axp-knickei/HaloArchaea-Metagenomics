# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2025-12-08

### Added
- **Initial Release**: Launched the HaloArchaea Metagenomics Pipeline, a high-throughput Snakemake v8 workflow designed for hypersaline ecosystems.
- **Workflow Orchestration**: 
    - Implemented support for Snakemake v8 with the `snakemake-executor-plugin-slurm` for cluster execution.
    - Added `main.py` CLI wrapper for simplified local and SLURM job submission.
- **Quality Control**: 
    - Integrated `fastp` with tuned parameters (relaxed complexity threshold) to preserve gas vesicle gene clusters (*gvp*) and other repetitive elements common in Halobacteria.
    - Added `MultiQC` aggregation for run reports.
- **Assembly**: 
    - Added `MetaSPAdes` (default) with multi-kmer iterative assembly for resolving strain heterogeneity.
    - Added support for `MEGAHIT` as an alternative assembler.
- **Binning**: 
    - Implemented `SemiBin2` with GPU acceleration support for deep learning-based binning.
    - Added `Bowtie2` mapping and sorting steps.
- **Taxonomy & Quality**: 
    - Integrated `GTDB-Tk` (v2.5.0) for taxonomic classification, featuring `skani` for rapid ANI screening.
    - Integrated `CheckM2` for machine learning-based genome quality assessment.
- **Annotation**: 
    - Added `DRAM` (Distilled and Refined Annotation of Metabolism) for functional annotation.
- **Specialized Analysis**: 
    - Added custom Python scripts (`calculate_pI.py`) to calculate proteome-wide isoelectric points (pI), validating the "salt-in" strategy typical of Haloarchaea.

### Configuration
- Defined resource profiles in `config/resources.yaml` for handling high-memory tasks (MetaSPAdes, DRAM) and GPU tasks (SemiBin2).
- Established container-first strategy using Apptainer/Singularity for major tools.