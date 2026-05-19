# OCM Demultiplexing - Development Journal

## 2026-05-06/07: Initial implementation (dev-ocm branch)

First OCM implementation on dev-ocm branch. Used pool-level demux CSV approach with large OB barcode files (921K entries each). Downstream QC ran per-pool via the custom demux path.

**Problem identified**: This approach incorrectly assumed pool-level ambient/doublet contamination. OCM (like Flex) physically separates samples, so ambient and doublets should be handled per-sample.

## 2026-05-08: Correct implementation (dev-flex branch)

Reimplemented OCM on dev-flex branch following the Flex pattern: map per pool → split to per-sample h5 → per-sample ambient/QC/doublets.

### Architecture refactor

- Merged `flex.smk` into `mapping.smk` with explicit `if IS_FLEX`/`else` branching
- `run_mapping` refactored: `{run}` → `{lib}` wildcard, params from `LIB_PARAMS`
- Unified `save_alevin_to_h5` handles non-muxed, Flex, and OCM (shared rule)
- `RUNS_TO_LIBS` maps sample → pool for both Flex and OCM

### Overhang-based barcode filtering

Discovered that OCM partition assignment is determined by a 2bp overhang at barcode position 8-9 (1-indexed):
- GT → OB1, CA → OB2, TC → OB3, AG → OB4
- Confirmed 100% against all 3.7M barcodes in 10x OB files
- Source: `overhang.txt` from cellranger (`lib/python/cellranger/barcodes/translation/overhang.txt`)

This eliminated the need for:
- Large OB barcode files (921K entries each)
- Manual downloads from 10x
- Chemistry-dependent file selection (3' vs 5')

Instead: `ocm_overhang_map.txt` is auto-extracted from cellranger during setup (~4 lines).

### Bug fixes during testing

- `tenx_assay_type` missing default in schema → added `"default": "poly_a"`
- `alevin_fry_home` path missing `ref_txomes/` subdirectory
- `mapping.py`: `tenx_chemistry_chem` variable name typo (UnboundLocalError)
- `mapping.py`: passed `tenx_chemistry` instead of `af_chemistry` to `_run_simpleaf_quant`
- `scprocess` CLI: `_render_index` called even after snakemake failure → added `pipefail` + exit check
- `mapping.py`: tmp fastq download dir now uses `lib_pool_dir` prefix
- `setup.py`: renamed `get_cellranger_whitelists` → `extract_cellranger_resources`

### Files modified

- `rules/mapping.smk` — merged flex.smk, unified save_alevin_to_h5, {lib} wildcard
- `rules/flex.smk` — gutted (content moved to mapping.smk)
- `rules/scprocess.smk` — IS_OCM, lib_pool_dir, unified af_mapping_outs
- `scripts/mapping.R` — .get_ocm_overhang helper, OCM filtering via substr
- `scripts/mapping.py` — lib_pool_dir support, bug fixes
- `scripts/scprocess_utils.py` — OCM validation, RUN_VAR/RUNS_TO_LIBS, ocm_id in params
- `scripts/SampleQC.R` — added "ocm" to demux_type checks
- `scripts/integration.py` — added "ocm" to demux_type checks
- `scripts/hvgs.py` — added "ocm" to demux_type checks
- `scripts/setup.py` — extract ocm_overhang_map.txt, renamed function
- `resources/schemas/config.schema.json` — "ocm" in enum, tenx_assay_type default
- `resources/mkdocs/docs/reference.md` — OCM docs
- `resources/mkdocs/docs/usage.md` — OCM usage section
- `scprocess` CLI — pipefail fix for render_index

## 2026-05-18: First OCM test run fixes

### save_alevin_to_h5 fixes (mapping.R)

The flex/OCM path in `save_alevin_h5_ambient_params` wrote the per-sample h5 and
set `fry_dir = NULL`, but then `save_alevin_h5_knee_params_df` tried to call
`loadFry(NULL)` → crash. Fixed by passing the in-memory matrix directly via
`precomputed_mat` parameter, avoiding h5 write→read round-trip entirely.

### Setup fixes

- Updated `TXOME_INDEX_FS` in `setup.smk`: prebuilt Zenodo tarballs now use
  `piscem_idx.ssi` (not `.sshash`) and `piscem_idx_ver.json` (not `piscem_idx.json`).
  The filename mismatch caused MissingOutputException, triggering snakemake to
  delete the successfully extracted files.
- Added download/extraction error checking in `_download_prebuilt_index`.

## 2026-05-19: Branch consolidation and bug fixes

### Branch merges

Merged `claude-shiny` (shiny app, join functionality) and `dev-fixes` (accumulated
bug fixes from main) into `dev`. Resolved conflicts where dev-flex refactoring
(run→lib, af_chemistry naming) diverged from the older branch structures.

### Chemistry variable naming cleanup

Renamed `sample_chem` → `sample_tenx_chemistry` in `mapping.py` for clarity.
Three chemistry variables now have distinct names:
- `tenx_chemistry` — user config value (10x kit name)
- `af_chemistry` — simpleaf `--chemistry` flag
- `sample_tenx_chemistry` — resolved per-sample chemistry at mapping time

### Zoom bug fixes

- `zoom.smk`: fixed `empty_gs_fs` passed as positional arg instead of `--empty_gs_f` flag to `hvgs.py calculate_hvgs`
- `zoom.smk`: fixed `ref_txome` → `genome_ref` parameter name in `run_fgsea` R call

### Documentation updates

- Added OCM setup note to `setup.md` (overhang map auto-extracted from cellranger)
- Expanded Flex and OCM descriptions in `usage.md`
- Fixed OCM reference in `reference.md` (overhang mechanism, not barcode files)
- Added `join.smk` to CLAUDE.md rule table
- Added chemistry variable naming convention to CLAUDE.md
