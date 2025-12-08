# Pipeline Details & Scientific Rationale

This pipeline is specialized for the metagenomic analysis of Class *Halobacteria* (Haloarchaea).

## Key Adaptations

### 1. "Salt-In" Strategy Verification
Haloarchaea accumulate molar concentrations of KCl to balance osmotic pressure. This results in a highly acidic proteome.
* **Pipeline Step:** `halo_analysis`
* **Method:** We extract all predicted proteins and calculate their isoelectric point (pI).
* **Expected Signature:** A median pI of ~4.0–5.0 confirms the organism uses the salt-in strategy.

### 2. Gas Vesicle Preservation
Haloarchaea often possess gas vesicle gene clusters (*gvp*) to regulate buoyancy. These regions are highly repetitive and often discarded by standard QC.
* **Tool:** `fastp`
* **Adaptation:** We relax the complexity threshold to **20** (default 30) to prevent over-trimming of these low-complexity regions.

### 3. Strain Heterogeneity
Hypersaline environments often contain multiple closely related strains.
* **Tool:** `MetaSPAdes`
* **Adaptation:** We use an iterative k-mer strategy (`k: [21, 33, 55, 77, 99, 127]`) to maximize graph connectivity and resolve microdiversity.

## Workflow Overview

| Stage | Tool | Notes |
| :--- | :--- | :--- |
| **QC** | fastp | Adapter trimming, poly-G tail trimming, complexity filtering. |
| **Assembly** | MetaSPAdes | De novo assembly. Uses high-memory nodes. |
| **Binning** | SemiBin2 | Deep learning-based binning. Uses GPU if available. |
| **Taxonomy** | GTDB-Tk | Classification against GTDB R220. Uses `skani` for fast ANI. |
| **Quality** | CheckM2 | Completeness/Contamination estimates. |
| **Annotation**| DRAM | Functional annotation and metabolic heatmap generation. |