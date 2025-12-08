# Contributing to HaloArchaea Metagenomics Pipeline

First off, thank you for considering contributing to the HaloArchaea Metagenomics Pipeline\! It's people like you that make this tool better for the scientific community.

We welcome contributions of all forms: bug reports, feature requests, documentation improvements, and code changes.

## 1\. How to Report Bugs

If you encounter a bug, please submit an issue on GitHub. To help us solve it quickly, please include:

  * **A clear title** describing the issue.
  * **Steps to reproduce** the bug.
  * **The command** you ran (e.g., `python main.py --mode slurm ...`).
  * **Log output** or error messages (check `logs/` directory).
  * Your system information (OS, Snakemake version).

## 2\. Setting Up Your Development Environment

To contribute code, you'll need to set up the environment locally.

### Prerequisites

  * [Mamba](https://github.com/mamba-org/mamba) or Conda.
  * Apptainer (Singularity) if you plan to test containerized rules.

### Installation

1.  **Fork** the repository on GitHub.

2.  **Clone** your fork locally:

    ```bash
    git clone https://github.com/YOUR-USERNAME/HaloArchaea-Metagenomics.git
    cd HaloArchaea-Metagenomics
    ```

3.  **Create the environment**:
    Follow the instructions in the `README.md` to create the master environment.

    ```bash
    mamba create -n snakemake_8 -c conda-forge -c bioconda \
        snakemake=8.4.6 \
        snakemake-executor-plugin-slurm \
        python=3.12 \
        pandas
    mamba activate snakemake_8
    ```

    *Note: Ensure you are using Python \>=3.12 as specified in `pyproject.toml`.*

4.  **Install the package**:

    ```bash
    pip install -e .
    ```

## 3\. Project Structure

Understanding the file layout will help you navigate the code:

  * **`workflow/Snakefile`**: The main entry point for the pipeline.
  * **`workflow/rules/`**: Modular Snakemake rules (e.g., `01_qc.smk`, `02_assembly.smk`).
  * **`workflow/envs/`**: Conda environment definitions for specific rules.
  * **`workflow/scripts/`**: Python scripts used by rules (e.g., `calculate_pI.py`).
  * **`config/`**: Configuration files (`config.yaml`, `resources.yaml`, `samples.tsv`).
  * **`main.py`**: The CLI wrapper for executing the workflow.

## 4\. Development Guidelines

### Snakemake Rules

  * **Modularity**: Keep rules separated by function in `workflow/rules/`.
  * **Resources**: Always define `resources` (mem\_mb, runtime) for computationally intensive rules. Use `config/resources.yaml` for defaults.
  * **Containers**: We prioritize a "Container-First" strategy. Ensure new rules specify a `container` directive.
  * **Logs**: Always include a `log` directive for rules that execute shell commands or scripts.

### Python Scripts

  * Scripts in `workflow/scripts/` should be executable and handle command-line arguments.
  * We recommend using `logging` instead of `print` statements, as seen in `calculate_pI.py`.
  * Follow PEP 8 style guidelines.

### Configuration

  * If you add a new tool, update `config/config.yaml` with its parameters.
  * If the tool requires significant resources (RAM/GPU), add a profile to `config/resources.yaml`.

## 5\. Testing

Before submitting a Pull Request, verify your changes.

1.  **Dry Run**: Ensure the DAG is correct.
    ```bash
    python main.py --dry-run
    ```
2.  **Local Test**: Run a small test on your local machine using the sample data provided (or dummy data).
    ```bash
    python main.py --mode local --cores 4
    ```

## 6\. Submitting a Pull Request

1.  Create a new branch for your feature or fix: `git checkout -b feature/amazing-feature`.
2.  Commit your changes with clear messages.
3.  Push to your fork: `git push origin feature/amazing-feature`.
4.  Open a Pull Request on the main repository.
5.  Describe your changes and link to any relevant issues.

## 7\. Documentation

  * If you change the pipeline usage, please update `README.md`.
  * If you change dependencies, update `pyproject.toml` and the relevant `workflow/envs/*.yaml` files.
