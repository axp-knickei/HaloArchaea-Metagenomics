# Database Setup

This pipeline relies on several large reference databases. You must download these before running the pipeline. Ensure you have ~200GB of disk space ready.

## 1. GTDB-Tk (Release 220)
Taxonomic classification relies on the Genome Taxonomy Database.

```bash
# Create directory
mkdir -p refs/gtdb_r220

# Download
wget [https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_r220_data.tar.gz](https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_r220_data.tar.gz)

# Extract
tar -xvzf gtdbtk_r220_data.tar.gz -C refs/gtdb_r220
````

*Update `config/config.yaml`*: Set `gtdb: "refs/gtdb_r220/release220"`

## 2\. CheckM2

Quality assessment uses a Diamond database.

```bash
# CheckM2 can download its own database
checkm2 database --download --path refs/checkm2_db
```

*Update `config/config.yaml`*: Set `checkm2: "refs/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd"`

## 3\. DRAM (Distilled and Refined Annotation of Metabolism)

DRAM requires multiple databases (KOFam, UniRef90, PFAM, etc.). This step requires high RAM (\>100GB).

```bash
DRAM-setup.py prepare_databases --output_dir refs/DRAM_data
```

*Update `config/config.yaml`*: Set `dram: "refs/DRAM_data"`
