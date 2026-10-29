# AGENTS.md

Instructions for coding agents working in this repository.

## Versioning, task scope and analysis provenance

### Required approval before main or release publication

Commits and pushes to `dev` are allowed within the authorized task. They do not
authorize promotion to `main` or a release. Before pushing/merging to `main`,
publishing any GitHub/Conda release or marking a version as a stable release,
obtain explicit approval from **both Will Macnair and Marusa** for the concrete
candidate commit, version and release contents. If either approval is absent or
unclear, leave the reviewed work on `dev` and ask; never infer Marusa's approval
from Will's approval, a generic request to fix/commit work, or passing tests.
Preparing and incrementing VERSION on `dev` is not release approval. Review
results must be available before seeking the two approvals. Material changes
after approval require renewed approval of the changed candidate.

Every completed version, including patch versions, must have an immutable
annotated Git tag `vX.Y.Z` pointing to the committed snapshot whose VERSION
matches. Will explicitly authorized development-version tagging on 2026-10-28.
These tags may be created and pushed on dev without publishing a release;
version tagging alone is not main promotion or stable-release approval.
Main promotion and release publication still require both maintainers' explicit
approval. Tag after validation and commit, never tag a dirty working tree as if
its uncommitted contents were included, and never move an existing shared tag.
The explicitly withdrawn v0.2.0/v0.2.1 tags remain absent.

Automatic publication on tag pushes is disabled. Do not re-enable publication
automation without an approval mechanism that requires both maintainers. The
user's explicit instruction to withdraw the premature releases and restore main
is authorization for that rollback only, not a new release.

Follow [the versioning and development policy](resources/mkdocs/docs/development/versioning.md).
At the start of a new task, inspect the checkout and active jobs, and establish
with the user whether software work is a patch, small additional functionality,
or major rewriting/breaking changes. Reuse a classification already agreed in
the conversation. For an analysis using unchanged software, identify the pinned
software version instead of inventing a release bump.

Keep implementation, meaningful tests and relevant docs in neat, compact commits.
Squash unpublished fixup/WIP commits when useful; never rewrite published,
shared, tagged or analysis-cited history without explicit authorization. Preserve
unrelated user edits. Prefer clean, committed checkouts for processing, training
and report runs; document any explicitly chosen dirty-run exception with a base
SHA, complete patch and relevant untracked source snapshots before execution.

`main` is the tested release branch; `dev` integrates development. For 0.x,
compatible fixes use patch releases, while features, breaks and intentional
scientific-default changes use minor releases. After 1.0, use normal semantic
versioning. VERSION, annotated release tag and GitHub release must agree. Record
breaking changes and changes affecting results separately. Models have immutable
independent IDs/checksums. Record analysis provenance per stage; a later commit
does not retroactively identify earlier dirty runs.


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
scprocess run <config.yaml> -r mapping|ambient|qc|hvg|integration|marker_genes|label_celltypes|zoom|shiny

# Unlock after an interrupted run
scprocess run <config.yaml> --unlock

# Create conda environments only (no execution)
scprocess run <config.yaml> --create-envs

# Build Shiny app from completed pipeline outputs
scprocess run <config.yaml> -r shiny

# Dry run the shiny build
scprocess run <config.yaml> -r shiny -n

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
| `label_celltypes.smk` | Cell type annotation (CellTypist or XGBoost); runs as part of `rule all` when configured |
| `train_xgboost.smk` | Train custom XGBoost classifier from annotated data; per-batch full-data evaluation |
| `render_htmls.smk` | HTML report generation from R Markdown |
| `zoom.smk` | Subclustering workflow for cell subsets |
| `shiny.smk` | Builds Shiny data, then performs a lightweight app deployment/configuration into `public/shiny_*`; invoked via `-r shiny` |
| `hto.smk` | Hashtag (multiplexing) demultiplexing; HTO index built with piscem via `simpleaf index` |
| `flex.smk` | Flex assay mapping (simpleaf multiplex-quant) |
| `join.smk` | Join multiple scprocess projects for combined analysis (BPCells PCA → Harmony) |
| `pb_empties.smk` | Pseudobulk ambient gene detection |

### Implementation Scripts (`scripts/`)

Python scripts are invoked by Snakemake rules using the `script:` directive. Key files:

- `scprocess_utils.py` — Config loading, JSON Schema validation, parameter resolution, setup verification, FASTQ metadata caching. Central to all rule setup.
- `hvgs.py` — Memory-efficient HVG detection: chunks the expression matrix and applies Seurat VST per chunk.
- `integration.py` — Scanpy-based integration: Harmony batch correction, optional GPU-accelerated RAPIDS, UMAP, Leiden clustering. Supports pre-computed PCA embeddings via `--precomputed_pca_f` (used by BPCells join path).
- `label_celltypes.py` — CellTypist and XGBoost-based cell type annotation. The `xgboost_one_batch` subcommand is also reused by `train_xgboost` for full-data evaluation.
- `train_xgboost.py` — Train XGBoost cell type classifier. Planning, fitting, aggregation, model identity helpers and label specification export live in this single Python file. Training-only imports are local so setup/config checks remain lightweight. Optional label_specs_f/coarse_label_specs_f CSVs control mapping, ordering and colours; their separate specs stage reuses saved predictions. The file starts with planning, then count preparation, fitting, output writing and evaluation. CLI entry functions retain public names; internal helpers use leading underscores. normalize_counts is shared with inference. Rules pass named validated settings rather than config_json. CLI uses subparsers: `prepare` (cell splitting and selection), `train` (model training) and `combine` (merge per-batch predictions with annotations and split labels) and `specs` (display tables and coarse evaluation). Works in scprocess, join, and zoom contexts — in join, h5ad lookup keys are `{project_id}_{batch_var}`.
- `mapping.py` — Chemistry auto-detection, simpleaf alignment, per-pool mapping. Main entry point is `map_fastqs_to_counts()`, which delegates to: `_resolve_fastq_paths()` (Arvados download or local path resolution), `_resolve_chemistry()` → `_autodetect_chemistry()` or `_validate_chemistry()` (chemistry/orientation detection), and `_write_chemistry_stats()` (saves results to YAML).

R scripts are invoked via `script:` directive or called in shell blocks:

- `mapping.R` — Converts alevin-fry output to per-sample h5 files; handles Flex (probe_id) and OCM (overhang) barcode splitting. For flex/OCM, the filtered matrix is passed in-memory via `precomputed_mat` to avoid h5 round-trip issues.
- `SampleQC.R` — QC metric visualization and filtering.
- `marker_genes.R` — Marker report plotting and expression helpers shared with classifier reports.
- `marker_genes_edger_bp.R` — The single edgeR.bp marker backend for projects, joins and zooms: BPCells pseudobulks, persisted TMM preparation, streamed fitting and compact report expression data.
- `label_celltypes.R` — R-side cell labeling utilities.
- `train_xgboost.R` — R plotting functions for the XGBoost training report: precision-recall plots, UMAP facets, gene importance dot plots, metrics tables.
- `shiny.R` — Builds and deploys the Shiny app. `make_shiny_app_scprocess()` reads scprocess-format inputs (`integrated_dt`, per-batch h5ads, markers, HVGs, fgsea) and writes BPCells-backed app data. `configure_shiny_app_scprocess()` performs lightweight annotation, template, and `shinyconfig.yaml` updates without rebuilding matrices.
- `join_pca.R` — BPCells disk-backed PCA for join workflow: converts the joint HDF5 counts once to persistent compressed `uint32` native BPCells storage, normalizes lazily (lib size → 1e4 → log1p), queues delayed per-gene centering/scaling, runs threaded `BPCells::svds()`, and writes PC scores to CSV.gz.

### Configuration System

Two-level config hierarchy:

1. **Setup config** (`scprocess_setup.yaml`, user-level): HPC profile, reference genomes, environment paths.
2. **Project config** (`config-{project}.yaml`): Sample metadata, QC thresholds, ambient removal method, integration settings, subclustering definitions.

Both configs are validated against JSON Schemas in `resources/schemas/` via `jsonschema` before any Snakemake execution. Schema files:
- `resources/schemas/config.schema.json` — project config
- `resources/schemas/setup.schema.json` — setup config
- `resources/schemas/zoom.schema.json` — zoom (subclustering) config
- `resources/schemas/join.schema.json` — join config
- `resources/schemas/custom_sample_params.schema.json` — per-sample overrides

### Execution Profiles

Snakemake profiles in `profiles/` define HPC-specific executor settings (job submission, GPU bindings, resource limits). Profiles: `slurm_default`, `lsf_default`, `slurm_shpc`.

Resource allocation for individual rules is predicted by linear models fitted on benchmark data from ~96 projects; parameters are stored in `resources/snakemake/resources_lm_params_2026-07-06.csv`. The prediction formula is `max(floor, intercept + slope × x)` where `x` is a rule-specific predictor (input file size, number of batches, etc.). Five per-sample rules (`save_alevin_to_h5`, `run_cellbender`, `get_barcode_qc_metrics`, `make_clean_h5ads`, `run_qc_one_run`) have separate model parameters for 10x v3 vs v4 chemistry. For auto-detected chemistry, the resolved chemistry is read from `chemistry_statistics.yaml` at resource-evaluation time. The CSV path is centralized via the `LM_PARAMS_FILENAME` constant in `scprocess_utils.py`.

### Conda Environments

Each rule specifies a `conda:` environment file from `envs/`. Key environments:
- `envs/scprocess_local.yaml` — main environment (Python, scanpy, polars, h5py, etc.)
- `envs/alevin_fry.yaml` — alignment
- `envs/rlibs.yaml` — R libraries (edgeR, fgsea, ggplot2)
- `envs/integration.yaml` — GPU-compatible scanpy/rapids stack
- `envs/bpcells_pca.yaml` — BPCells + irlba for disk-backed PCA (join workflow)
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
→ [scprocess run -r shiny] → BPCells matrix + shinyconfig.yaml → public/shiny/
```

The primary data format throughout is `AnnData` (`.h5ad`), with intermediate SCE objects (`.rds`) for R-based steps.

The `shiny` rule (`-r shiny`) requires a completed `marker_genes` run (including GSEA). It reads `integrated_dt_*.csv.gz` for UMAP coordinates and cluster labels (`RNA_snn_res.{mkr_sel_res}` column), combines per-batch h5ad files via zellkonverter for the count matrix, and writes a density-subsampled BPCells matrix to `public/shiny/data/`.

The project config may include an optional `shiny:` section to customise `app_title`, `email`, `keyword`, `default_gene`, `n_keep` (cells to subsample, default 30000), `metadata_vars` (override which variables to show), `metadata_labels` (display names as key-value pairs), and `metadata_combns` (variable pairs to cross-tabulate).

## Chemistry Detection

The pipeline determines 10x chemistry at two stages:

### Config-time (`scprocess_utils.py:_get_lib_parameters_one_lib()`)

`tenx_chemistry` is set in the project config (default `"none"` = auto-detect). Maps to `af_chemistry` (the alevin-fry `--chemistry` flag):

- `3v2`, `5v1`, `5v2` → `af_chemistry = '10xv2'`, shared 737K whitelist
- `3v3`, `3v4`, `3LT`, `5v3`, `multiome` → `af_chemistry = '10xv3'`, chemistry-specific whitelist
- `none` → auto-detect at mapping time

Orientation: 3' chemistries → `fw`; 5' chemistries → `rc`. Whitelist files looked up from `$SCPROCESS_DATA_DIR/cellranger_ref/cellranger_whitelists.csv`. Per-sample overrides possible via `custom_sample_params.yaml`.

### Mapping-time auto-detection (`scripts/mapping.py:_autodetect_chemistry()`)

When `tenx_chemistry == "none"`, three-step auto-detection:

1. **Barcode overlap** (`_get_whitelist_overlap()`) — samples 100K R1 reads, extracts 16bp barcodes, tests overlap against all whitelists. Picks highest overlap (warns if < 70%).
2. **R1 length check** — if `10xv3` but R1 < 28bp, falls back to `10xv2` with v3 whitelist (warns).
3. **Orientation inference** — if multiple chemistries share a whitelist (e.g. `3v2`/`5v1`/`5v2`), downsamples FASTQs, maps twice (fw + rc), picks orientation with more quantified cells.

When `tenx_chemistry` is set explicitly, `_validate_chemistry()` checks R1 length compatibility and falls back to `10xv2` with a warning if R1 < 28bp.

Results saved per-run to `chemistry_statistics.yaml` via `_write_chemistry_stats()`, merged into `chemistry_statistics_all_runs_{DATE_STAMP}.csv`.

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

### Mapping output folder structure

Subdirectories under `af_{sample_or_pool}/` depend on the assay type:

- **Simple (no multiplexing)**: `af_{sample}/af_quant/`, `af_{sample}/af_counts_mat.h5` — no subdirectory
- **HTO**: `af_{pool}/rna/af_quant/` and `af_{pool}/hto/af_quant/` — `rna/` subdirectory only exists to separate RNA from HTO outputs
- **Flex**: `pools/af_{pool}/flex/af_quant/` — under `pools/` prefix with `flex/` subdirectory

Controlled by `af_rna_dir` in `scprocess.smk`: `'flex/'` for flex, `'rna/'` for HTO, `''` otherwise.

## Coding Conventions

The shared user preferences are in
[/home/macnairw/.codex/coding-style.md](/home/macnairw/.codex/coding-style.md).
Read that guide when available; the repository-specific conventions below
apply to scprocess. Review history and unfinished batches remain in
`../scprocess-design/style_review_plan.md` and `style_review_examples.md`.

### Implementation layout

- Use established two-space Python/Snakemake indentation and local assignment/
  argument alignment, including spaces around keyword-argument equals signs.
  Do not impose one quote style across the repository.
- Keep the established R idioms: function definitions use `<-`, ordinary local
  assignments generally use `=` with local alignment, and data.table and `%>%`
  remain appropriate. Do not change unusual assignment contexts mechanically.
- Keep related small training helpers in train_xgboost.py rather than splitting
  off tiny planning/workflow/model modules; normalize_counts remains shared with
  inference. This is the reviewed organisation of the current trainer.
- When a script runs in multiple Conda environments, keep genuinely heavy or
  training-only imports local where needed, with a brief reason. This does not
  justify repeated local NumPy imports throughout one file.
- Pass explicit project/join context (is_join) to classifier planning; batch_var
  is a column name and may be sample_id or pool_id. Use integration_f consistently
  for the integration table in scripts, CLI flags and named workflow inputs.

### Schema defaults and CLI settings

User clarification, 2026-10-13: never introduce fallback defaults for configurable
parameters in workflow rules or implementation helpers. In particular, prohibit
config.get('parameter', literal_default), config['section'].get('parameter',
literal_default), setdefault with a literal default, and equivalent expressions
that substitute a value when a validated setting is missing. Put default values
in the relevant schema and apply them during configuration validation. Consumers
must access the validated settings explicitly. A missing required value is a
validation error, not permission to choose another default locally.

Optional-feature checks may represent an inactive argument as absent/NULL;
these must not choose an alternative active-analysis value. Conditional defaults
must remain conditional, with values defined in the schema. Generic lookups for
optional sections and computed file paths are not parameter default definitions.
Review existing local fallback defaults as findings in their own batches rather
than silently expanding a style edit into a repository-wide behavior change.

The schemas are in resources/schemas/. Essential file paths, identifiers and
schema-supplied scientific settings remain explicit CLI inputs. Optional features
may be omitted; check paired options together. CLI declaration groups should be
easy to scan by type. Relative paths resolve from the project folder; data/configs/
is a location used in examples, not an implicit filename search directory.

Fixed method constants are distinct from configurable defaults. For example,
MIN_CLUSTER_MAJORITY_VOTE_CELLS is the fixed 10-training-cell minimum for a cluster
vote, not a YAML option. Keep such constants named and explained at file level.
Classifier aggregation requires the training-metadata JSON to determine whether legacy
coarse mapping was requested; missing metadata is an error, while a recorded
absence of optional mapping is valid.

### Snakemake workflow files

Treat workflow files as a distinct review area. `rules/join_zoom.smk` is an
historical negative example identified by the user, not a style reference.
Compare it with established setup and rule patterns in historical integration,
zoom, HVG and marker workflows. Review:

- Setup sequence: imports, config validation, parameter extraction, common
  variables/output paths, targets and rules should be easy to follow.
- Names and scope: distinguish parameter dictionaries, lists of batches, paths
  and constants. Avoid opaque aliases such as `_z` and unexplained context
  globals such as `JOIN_ZOOM_CONTEXT`; document and reduce implicit coupling
  across `include` boundaries where a concrete redesign supports that.
- Use uppercase names for file-level code constants such as
  `JOIN_ZOOM_SCRIPT_F` and `HVGS_SCRIPT_F`. All named rule inputs use descriptive
  lowercase names, for example `join_zoom_script_f = JOIN_ZOOM_SCRIPT_F`.
  Code remains a tracked input.
- Give rules descriptive names that identify their calculation and output;
  update corresponding resource keys, profile references and tests together.
- Declare file dependencies under `input:`, generated files under `output:`,
  and individual named settings under `params:`. Lists and dictionaries may be
  passed directly to `script:` rules. Avoid generic JSON `context`/`effective`
  bundles. Use the supplied join/zoom YAML files as the configuration record;
  do not add a separate resolved-configuration snapshot.
- Keep rules readable and primarily declarative. Put substantial validation,
  selection and serialization logic in appropriately named helpers; do not
  hide simple path declarations behind unnecessary helper layers.
- Use established two-space indentation, local alignment, readable multiline
  directives and shell arguments, and short section comments. Compare directive
  order with adjacent historical rules rather than blindly reformatting all files.
- Rule names and ordering should communicate the analytical stages and their
  dependencies.
- Preserve logging, resources, environment choices, wildcard behavior and DAG
  dependencies. Keep `script:` where it solves a real transport/execution need;
  do not replace it with `shell:` merely for visual consistency.
- Any restructuring must preserve rerun behavior, tracked configuration, source
  isolation, optional outputs and the distinction between join and joined zoom.
- Keep .smk files flat and declarative. Configuration checks live in the
  _check_*_parameters() functions in scprocess_utils.py, not in rule files.
- Prefer unconditional rule definitions with optional targets controlled by the
  DAG, using empty lists where appropriate rather than if blocks around rules.
- Follow directive order input → output → params → threads → retries → resources
  → log → benchmark → conda → shell, preserving justified adjacent patterns.
- Shell log redirection starts with exec &>> {log}. Preserve scoped redirection
  for local Python execution and justified script: interfaces.
- Do not pass generic config_json bundles. Changing unrelated configuration
  must not force expensive upstream recomputation.

### Reports and documentation

Review reports independently from the website/reference documentation. The
explanatory text is for computational biologists. For each substantive output
section, assess whether it:

1. Provides biological context: what question the figure/table addresses and
   what its axes, groups, colours or summary represent.
2. Indicates success and failure: recognizable evidence of a useful result,
   warning signs, and limitations of interpretation. Distinguish a technical
   failure from a biologically plausible absence of signal; do not prescribe
   an attractive plot as proof of a valid analysis.
3. Briefly explains parameter adjustment: name the relevant configuration
   parameter, explain what it affects and when changing it may help. Verify
   names and behavior against the actual schema/code. If no adjustment can
   address the concern, say that rather than inventing a tuning suggestion.

These are content criteria, not a requirement to add three repetitive headings
to every plot. Prefer concise, connected explanations near the outputs.
Avoid explaining worker pools, serialization, caching mechanics, implementation
history, seeded display draws or source-code organization unless they change
how a computational biologist should interpret or use the result. Retain
necessary methodological details such as normalization, transformations and
scaling when they change the meaning of an axis or comparison. Compact method
notes may support the explanation; developer details belong in developer notes.

Review figure captions, tab introductions, missing-result messages and prose
generated by R chunks as well as static Markdown. Include all report stages:
mapping, ambient removal, multiplexing, QC, HVGs, integration, cell cycle,
markers, annotation/training, join, zoom and index summaries. Evaluate each
according to its purpose rather than inserting boilerplate everywhere.

Reuse established integration.Rmd.template explanations for corresponding join
outputs, checking historical and related templates before drafting new text.
Match text to the actual plot: entropy guidance does not describe a composition
UMAP. Record the reused source section and any real gap during a report review.
The superseded agent-generated join-integration preview in the examples document
is historical evidence, not an approved style reference.

Keep reusable preparation and plotting functions in associated R helpers and
heavy calculations outside rendering. Validate representative rendered HTML after
report changes. Website, schema, CLI and Shiny text have their own audiences:
retain actionable setup, resource and dependency guidance where it belongs.
Check resources/mkdocs when functionality changes. Preserve the user's existing
cell-cycle introduction and unrelated edits during review.

## Key Patterns

- **Config access in scripts**: All scripts receive `snakemake.config` and `snakemake.params`; utilities are imported from `scprocess_utils.py` for consistent parameter resolution.
- **Chunked processing**: HVG detection and some QC steps process data in chunks to handle datasets with millions of cells.
- **Batch iteration**: Rules in `scprocess.smk` expand over `batches` defined in the project config, enabling parallel sample processing.
- **HTML reports**: Rendered via `render_htmls.smk` from R Markdown templates in `resources/rmd_templates/`; generated at QC, integration, and marker gene stages.
- **BPCells PCA for join**: The join workflow supports disk-backed PCA via BPCells (`int_pca_method: bpcells`, default) to handle datasets exceeding GPU memory. A separate rule converts the joint HDF5 count matrix to persistent compressed `uint32` native BPCells storage, which is copied to local `/tmp` for fast I/O. Per-gene centering/scaling stay lazy through BPCells operator transformations, and threaded `BPCells::svds()` performs the iterative SVD. The resulting PC embeddings are passed to Harmony on GPU.
- **Join rule separation**: `join_build_matrix` (expensive, loads h5ads) is separate from `join_build_coldata` (lightweight, reads barcodes from HDF5). This means changing `metadata_vars` only reruns coldata + render, not the full matrix assembly. PCA reads library sizes directly from the matrix via `BPCells::colSums()`, independent of coldata.
- **label_celltypes in rule all**: When a `label_celltypes` block is present in the config, its outputs are included in `rule all` automatically. The `save_cluster_names_file: true` option generates a `cluster_names_for_shiny_*.csv` at the `marker_genes:mkr_sel_res` resolution.
- **QC checkpoint for label_celltypes**: `get_qc_sample_statistics` is a Snakemake `checkpoint` (not a regular rule). The `merge_labels` rule uses an input function (`_get_good_batch_labels`) that reads the checkpoint output to discover which batches passed QC. This allows label_celltypes to run in a single `scprocess run` invocation without requiring QC to have completed before DAG construction.
- **FASTQ metadata caching**: `_get_fastqs()` writes per-library YAML cache files (`fastq_metadata_{lib}.yaml`) to the mapping output directory after resolving FASTQ files (from filesystem or Arvados). Subsequent runs load from cache, allowing FASTQ files to be deleted after mapping without breaking downstream rules. The `_check_fastq_dir_or_cache()` helper tolerates a missing `fastq_dir` when cached metadata exists.
- **R environment isolation**: `_clean_conda_env()` in the CLI sets `R_LIBS_SITE=""`, `R_LIBS_USER=""`, `BASILISK_USE_SYSTEM_DIR=0`, and `BASILISK_EXTERNAL_DIR` to prevent conda R from loading system packages and to store basilisk Python envs alongside Snakemake conda envs. R scripts that use `zellkonverter` via parallel workers (e.g. `marker_genes.R`) pre-initialize the basilisk env before `bplapply` to avoid race conditions.
- **Join label reuse**: The join workflow checks whether source projects already have label_celltypes outputs for the same labeller/model. If they do, naive predictions are extracted from the source labels instead of re-running celltypist/XGBoost — only projects without existing labels trigger fresh runs.
- **Shiny as a rule**: Shiny is invoked via `scprocess run -r shiny` or `scprocess join -r shiny` (not a separate subcommand). `shiny.smk` separates expensive `build_*_shiny_data` rules from lightweight `configure_*_shiny_app` rules; the CLI targets configuration, which depends on data generation only when analytical inputs are stale.
- **train_xgboost across contexts**: `train_xgboost` runs in three contexts (scprocess, join, zoom), each with different h5ad key construction. The preparation rules pass explicit is_join context: false for project/project zoom and true for join/joined zoom. Do not infer this from batch_var. Project H5AD keys match directly; joined keys are project-qualified, preserving already-qualified integration values. All shell blocks must use the `train` subcommand (`python3 scripts/train_xgboost.py train ...`). Full-data evaluation reuses `label_celltypes.py xgboost_one_batch` for per-batch prediction, with `combine` merging results and assigning train/validation/holdout/unlabeled split labels.
