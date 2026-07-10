# UCI Single-Cell RNA-seq Pipeline

This repository provides a comprehensive Snakemake-based workflow for analyzing single-cell RNA-seq data, optimized for 10x Genomics datasets and HPC environments using SLURM.

## Main Features
- **CellRanger Processing:** Automated workflows for single and multi-sample analysis using CellRanger.
- **Downstream Analysis:** Includes rules and scripts for SCENIC (gene regulatory network), Numbat (CNV analysis), SNP phasing, and genome browser visualization.
- **Modular Design:** Workflows are organized into rule files and scripts for flexible, reproducible analysis.
- **HPC Support:** Resource requirements and SLURM integration for large-scale data processing.

![The workflow](rulegraph.png)

## Directory Structure
- `Snakefile_cellranger`, `Snakefile_cellranger_multi`: Main Snakemake workflows.
- `config_cellranger.yaml`, `config_cellranger_multi.yaml`: Configuration files for sample lists and paths.
- `rules/`: Snakemake rule modules for downstream analysis.
- `scripts/`: R, Python, and shell scripts for data processing and visualization.
- `output/`: Results and reports from pipeline runs.
- `data/`: Input FASTQ files and sample metadata.
- `envs/`: Conda environment YAMLs for reproducibility.

## Getting Started
0. **Clone the repository:**
   ```bash
   git clone https://github.com/uci-grthub/uci_grthub_parse_pipeline.git
   cd uci_grthub_parse_pipeline
   ```
1. **Install Snakemake:**
   ```bash
   conda install -c conda-forge snakemake
   ```
   The downstream Python/R analysis stack (scanpy, Harmony, Seurat, Shiny, etc.)
   is pinned separately in `pixi.toml`; install it with:
   ```bash
   pixi install -e uci-grthub-parse-pipeline
   ```

2. **Configure YAML files:** Edit sample IDs and paths in `config_cellranger.yaml` or `config_cellranger_multi.yaml`.
3. **Run Dry-Run Validation:**
   ```bash
   snakemake --snakefile Snakefile_cellranger --configfile config_cellranger.yaml --dryrun
   ```
4. **Execute Workflow:**
   - Local: `snakemake --snakefile Snakefile_cellranger --configfile config_cellranger.yaml --cores 8`
   - SLURM: `sbatch submit_cellranger.sh`

## Requirements
- CellRanger 8.0.1
- R/Bioconductor packages (install via `scripts/install_pkgs.R`)
- Reference genomes in `/dfs8/commondata/cellranger/`
- SLURM cluster environment

## Troubleshooting
- Ensure FASTQ files follow 10x naming conventions.
- Always validate configuration and paths before running.
- Check log and benchmark files for errors.

## Documentation
- See `README_cellranger.md` for CellRanger workflow details.
- See `rules/README.md` for downstream analysis rules.

For more information, contact the pipeline maintainers or consult the documentation files in this repository.
