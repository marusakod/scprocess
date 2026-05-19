# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project-specific Instructions for Development

- When you have finished working on some functionality, you MUST check that the docs (in resources/mkdocs) are up-to-date.

## Project Overview

**scprocess** is a Snakemake pipeline for automated processing of single-cell and single-nuclei RNA-seq data from 10x Genomics technology. It orchestrates a series of bioinformatics steps from raw reads through QC, ambient RNA removal, dimensionality reduction, batch correction, and marker gene detection.

Stack: Python 3.12 + R + Snakemake 9.8.1, with conda environments managed per rule.

## Common Commands

```bash
# Run the full pipeline for a project
scprocess run <config.yaml>

# Dry run (show what would execute without running)
scprocess run <config.yaml> -n

# Run a specific pipeline stage
scprocess run <config.yaml> -r mapping|ambient|qc|hvg|integration|marker_genes|label_celltypes|zoom

# Unlock after an interrupted run
scprocess run <config.yaml> --unlock

# Create conda environments only (no execution)
scprocess run <config.yaml> --create-envs

# Build Shiny app from completed pipeline outputs
scprocess shiny <config.yaml>

# Dry run the shiny build
scprocess shiny <config.yaml> -n

# Generate interactive knee plot for a sample
scprocess plotknee <sample> -c <config.yaml>

# Create new project directory scaffold
scprocess newproj <name> -w <where> [-s] [-c sc|sn|multiplex]

# Build docs site (in resources/mkdocs/)
mkdocs serve  # local preview
mkdocs build  # build static site
```

## Architecture

### Entry Point

`scprocess` (root-level executable) — a Python CLI dispatcher for `setup`, `newproj`, `run`, and `plotknee` subcommands. It constructs and invokes the appropriate `snakemake` command with the correct profile and flags.

### Workflow Rules (`rules/`)

Each `.smk` file is a self-contained Snakemake rule module. The core orchestrator is `rules/scprocess.smk`, which imports other modules and defines the DAG:

| Rule file | Responsibility |
|-----------|---------------|
| `scprocess.smk` | Core orchestrator; imports all modules |
| `setup.smk` | Reference genome indexing |
| `mapping.smk` | Read alignment via simpleaf (alevin-fry) |
| `ambient.smk` | Ambient RNA removal (CellBender or DecontX) |
| `qc.smk` | Cell QC, doublet detection (scDblFinder), filtering |
| `hvgs.smk` | Highly variable gene detection (chunk-wise VST) |
| `integration.smk` | Harmony batch correction, PCA/UMAP, clustering |
| `marker_genes.smk` | Pseudobulk DE with edgeR + GSEA |
| `label_celltypes.smk` | Cell type annotation (CellTypist or XGBoost) |
| `render_htmls.smk` | HTML report generation from R Markdown |
| `zoom.smk` | Subclustering workflow for cell subsets |
| `shiny.smk` | Standalone rule: builds and deploys Shiny app into `public/shiny/` |
| `hto.smk` | Hashtag (multiplexing) demultiplexing |
| `flex.smk` | Flex assay mapping (simpleaf multiplex-quant) |
| `join.smk` | Join multiple scprocess projects for combined analysis |
| `pb_empties.smk` | Pseudobulk ambient gene detection |

### Implementation Scripts (`scripts/`)

Python scripts are invoked by Snakemake rules using the `script:` directive. Key files:

- `scprocess_utils.py` — Config loading, JSON Schema validation, parameter resolution, setup verification. Central to all rule setup.
- `hvgs.py` — Memory-efficient HVG detection: chunks the expression matrix and applies Seurat VST per chunk.
- `integration.py` — Scanpy-based integration: Harmony batch correction, optional GPU-accelerated RAPIDS, UMAP, Leiden clustering.
- `label_celltypes.py` — CellTypist and XGBoost-based cell type annotation.
- `mapping.py` — Chemistry auto-detection, simpleaf alignment, per-pool mapping.

R scripts are invoked via `script:` directive or called in shell blocks:

- `mapping.R` — Converts alevin-fry output to per-sample h5 files; handles Flex (probe_id) and OCM (overhang) barcode splitting. For flex/OCM, the filtered matrix is passed in-memory via `precomputed_mat` to avoid h5 round-trip issues.
- `SampleQC.R` — QC metric visualization and filtering.
- `marker_genes.R` — edgeR pseudobulk DE testing and GSEA visualization (fgsea).
- `label_celltypes.R` — R-side cell labeling utilities.
- `shiny.R` — Builds and deploys the Shiny app; exports `make_shiny_app_scprocess()`. Reads scprocess-format inputs (`integrated_dt`, per-batch h5ads, markers, HVGs, fgsea), writes a BPCells on-disk count matrix, and produces `shinyconfig.yaml` in the deploy directory.

### Configuration System

Two-level config hierarchy:

1. **Setup config** (`scprocess_setup.yaml`, user-level): HPC profile, reference genomes, environment paths.
2. **Project config** (`config-{project}.yaml`): Sample metadata, QC thresholds, ambient removal method, integration settings, subclustering definitions.

Both configs are validated against JSON Schemas in `resources/schemas/` via `jsonschema` before any Snakemake execution. Schema files:
- `resources/schemas/config.schema.json` — project config
- `resources/schemas/setup.schema.json` — setup config
- `resources/schemas/zoom.schema.json` — zoom (subclustering) config
- `resources/schemas/custom_sample_params.schema.json` — per-sample overrides

### Execution Profiles

Snakemake profiles in `profiles/` define HPC-specific executor settings (job submission, GPU bindings, resource limits). Profiles: `slurm_default`, `lsf_default`, `slurm_shpc`.

Resource allocation for individual rules is predicted by an ML model; parameters are stored in `resources/snakemake/resources_lm_params_2025-12-16.csv`.

### Conda Environments

Each rule specifies a `conda:` environment file from `envs/`. Key environments:
- `envs/scprocess_local.yaml` — main environment (Python, scanpy, polars, h5py, etc.)
- `envs/alevin_fry.yaml` — alignment
- `envs/rlibs.yaml` — R libraries (edgeR, fgsea, ggplot2)
- `envs/integration.yaml` — GPU-compatible scanpy/rapids stack
- `envs/shiny.yaml` — Shiny app build environment (BPCells, zellkonverter, shiny, shinydashboard, ComplexHeatmap, Seurat, etc.)

### Shiny App (`resources/shiny/`)

Template files for the deployed Shiny app:
- `app.R`, `constants.R` — entry point and shared constants
- `utils/` — data I/O, colour palettes, plot helpers
- `modules/` — four tab modules: `explore_genes.R`, `explore_clusters.R`, `explore_genesets.R`, `explore_prevalence.R`
- `extdata/genesets/` — must be populated by the user with `transcription_factors_{human,mouse}.txt.gz` and `genes_go_pathways_{human,mouse}.txt.gz` (copied from the scjoin package)

The deployed app reads `shinyconfig.yaml` (written by `scripts/shiny.R`) to locate data files. Species is inferred from `ref_txome` (GRCh/hg → human; GRCm/mm → mouse).

## Data Flow

```
Raw FASTQ → simpleaf alignment → SCE (per sample)
→ Ambient RNA removal → Cell calling → Doublet detection
→ QC filtering → HVG detection (VST)
→ Harmony integration → PCA → UMAP → Leiden clustering
→ Marker genes (edgeR + GSEA) → HTML reports
→ [scprocess shiny] → BPCells matrix + shinyconfig.yaml → public/shiny/
```

The primary data format throughout is `AnnData` (`.h5ad`), with intermediate SCE objects (`.rds`) for R-based steps.

The `shiny` subcommand requires a completed `marker_genes` run (including GSEA). It reads `integrated_dt_*.csv.gz` for UMAP coordinates and cluster labels (`RNA_snn_res.{mkr_sel_res}` column), combines per-batch h5ad files via zellkonverter for the count matrix, and writes a density-subsampled BPCells matrix to `public/shiny/data/`.

The project config may include an optional `shiny:` section to customise `app_title`, `email`, `keyword`, `default_gene`, `n_keep` (cells to subsample, default 30000), `var_names`, and `var_combns`.

## Chemistry Detection

The pipeline determines 10x chemistry at two stages:

### Config-time (`scprocess_utils.py:_get_lib_parameters_one_lib()`)

`tenx_chemistry` is set in the project config (default `"none"` = auto-detect). Maps to `af_chemistry` (the alevin-fry `--chemistry` flag):

- `3v2`, `5v1`, `5v2` → `af_chemistry = '10xv2'`, shared 737K whitelist
- `3v3`, `3v4`, `3LT`, `5v3`, `multiome` → `af_chemistry = '10xv3'`, chemistry-specific whitelist
- `none` → auto-detect at mapping time

Orientation: 3' chemistries → `fw`; 5' chemistries → `rc`. Whitelist files looked up from `$SCPROCESS_DATA_DIR/cellranger_ref/cellranger_whitelists.csv`. Per-sample overrides possible via `custom_sample_params.yaml`.

### Mapping-time auto-detection (`scripts/mapping.py`)

When `tenx_chemistry == "none"`, three-step auto-detection:

1. **Barcode overlap** — samples 100K R1 reads, extracts 16bp barcodes, tests overlap against all whitelists. Picks highest overlap (warns if < 70%).
2. **R1 length check** — if `10xv3` but R1 < 28bp, falls back to `10xv2` with v3 whitelist.
3. **Orientation inference** — if multiple chemistries share a whitelist (e.g. `3v2`/`5v1`/`5v2`), downsamples FASTQs, maps twice (fw + rc), picks orientation with more quantified cells.

Results saved per-run to `chemistry_statistics.yaml`, merged into `chemistry_statistics_all_runs_{DATE_STAMP}.csv`.

### Chemistry variable naming

Three related variables are used across the codebase:

- `tenx_chemistry` — user-facing config value (10x kit name: `3v3`, `5v2`, `flexv1`, `none`, etc.)
- `af_chemistry` — simpleaf `--chemistry` flag (`10xv2`, `10xv3`, `10x-flexv1-gex-3p`, etc.)
- `sample_tenx_chemistry` — the resolved `tenx_chemistry` for a specific sample at mapping time, written to YAML as `selected_tenx_chemistry`

## Multiplexing Architecture

Three demultiplexing modes, each with a different data flow:

| Mode | Mapping unit | H5 unit | Demux mechanism |
|------|-------------|---------|-----------------|
| `none` | sample | sample | N/A |
| `hto` | pool | pool | Seurat HTODemux on HTO counts |
| `ocm` | pool | sample | Overhang-based barcode splitting (2bp at position 8-9) |
| `flex` | pool | sample | Probe ID prefix filtering |
| `custom` | pool | pool | User-supplied demux CSV |

For Flex and OCM, pool-level mapping output is split into per-sample h5 files at the `save_alevin_to_h5` rule, so ambient/QC/doublet processing runs per-sample. The `{lib}` wildcard handles pool-level mapping; `{run}` handles per-sample downstream processing.

## Key Patterns

- **Config access in scripts**: All scripts receive `snakemake.config` and `snakemake.params`; utilities are imported from `scprocess_utils.py` for consistent parameter resolution.
- **Chunked processing**: HVG detection and some QC steps process data in chunks to handle datasets with millions of cells.
- **Batch iteration**: Rules in `scprocess.smk` expand over `batches` defined in the project config, enabling parallel sample processing.
- **HTML reports**: Rendered via `render_htmls.smk` from R Markdown templates in `resources/rmd_templates/`; generated at QC, integration, and marker gene stages.
