# Changelog

## 0.3.0 — Unreleased development candidate

### Training features and interface changes

- Add saved XGBoost training plans, explicit test groups, whole-group validation,
  joint variable-combination sampling and separate training/validation caps.
- Rename training `ref_tag` to `xgboost_tag`; replace `refine_labels` with optional
  `cluster_majority_vote_resolution`. Omitted resolution preserves original
  labels, changing the former smoothing default. Rename the associated purity
  option and reject obsolete YAML keys; see the XGBoost training guide.
- Compute majority votes from training annotations only, preserve original test
  labels, filter genes using training cells only, align reordered genes by ID,
  and slice exported/predicted models at the best early-stopping iteration.
- Separate independent test results from early-stopping validation and sampled
  development holdouts in full-data summaries and reports.
- GPU training selects an explicit CUDA-enabled environment and verifies real GPU
  execution before reading expression; missing CUDA/CPU fallback is an error.
- Candidate only: no release, main promotion or model registration is implied.

## 0.2.2 — Unreleased development version

- Require explicit approval from both Will Macnair and Marusa before promotion
  to main, release tags or release publication. Dev commits/pushes and VERSION
  updates do not imply release approval. Disable automatic release publication.
- Retain recent processing and Shiny fixes on dev for review. Restore main to
  the previous stable processing code, retaining the approval policy and docs CI.
- The 0.2.0 and 0.2.1 publications below were premature and withdrawn at the
  user's request; preserve their commits and analysis provenance.
- Next: train/evaluate/register the Siletti classifier as a separate feature
  version; then resume project-specific code/docs/Rmd style and tone review.


## 0.2.1 — 2026-10-25 (publication withdrawn)

### Fixes

- Keep explicit read-only HDF5 handles alive through Shiny's lazy BPCells reads.
  This prevents BPCells 0.3.1 from modifying source H5AD timestamps. Add hdf5r
  explicitly to the managed Shiny environment. Counts and analysis parameters
  are unchanged; previously touched timestamps are not reset.
- Update GitHub Pages to Python 3.12 on Ubuntu 24.04. The previous Python 3.9
  setup stopped working after ubuntu-latest migrated to Ubuntu 26.04. This fixes
  documentation deployment without altering the immutable 0.2.0 release.

### Validation

- Regression coverage checks count preservation, lazy cross-file evaluation,
  source checksums and mtimes, and closure of read-only handles. Existing Shiny
  gene-alignment, density and workflow/configuration checks are retained.
- Siletti source inputs are checked against the pre-Shiny filesystem snapshot
  before refreshing failed-build input signatures; timestamp changes alone are
  not silently accepted as proof of unchanged contents.


## 0.2.0 — 2026-10-25 (publication withdrawn)

This release promotes the accumulated dev work since the published GitHub
v0.1.1 release to main. It establishes the versioning/development policy and
records the outstanding Siletti-related fixes. See the comparison with v0.1.1
for the complete change history; this is not a claim that all Siletti processing
stages used exactly this final release.

### Breaking changes and migration

- Doublet-cluster threshold methods are now `adaptive` and `manual`. Migrate
  `split_mad_mixture` to `adaptive` and `fixed` to `manual`. `mean_mad` and
  `int_dbl_mad_multiplier` were removed.
- `adaptive` rejects an explicit `int_dbl_cl_prop`; remove it from adaptive
  configs. `manual` requires an explicit `int_dbl_cl_prop`, with no default.
- Joined marker preparation/calculation rules and disk-backed PCA stages have
  been reorganized. Update scripts targeting internal historical rule names;
  use the documented stage interfaces. Review current schemas and release
  documentation before migrating saved configs or resuming old workflows.

### Changes affecting scientific results

- Adaptive filtering uses a robust split median/MAD logit mixture over preliminary
  cluster doublet fractions. This differs from historical fixed/mean-MAD rules;
  preserve the saved effective cutoff and excluded clusters for each run.
- scDblFinder uses package-default options unless explicitly configured; optional
  per-run uncertainty overrides, cluster-based calling and a separate UMI input
  filter are supported. These settings can change doublet calls.
- The supported scDblFinder/XGBoost environment is pinned. A call-local scorer
  compatibility fix aligns predictor rows after artificial-doublet filtering;
  training failures raise an error instead of accepting fallback scores.
- Cell-cycle estimation/correction and integration options are explicit. Review
  current PCA, Harmony, GPU and cell-cycle choices rather than assuming that an
  older configuration describes identical processing.

### Features

- Expanded join and joined-zoom workflows, BPCells PCA and streamed pseudobulk
  marker fitting, optional cell-cycle correction and multi-covariate integration.
- Upper QC count/feature limits, richer doublet diagnostics, mapping knee/shin
  scatter diagnostics and updated QC/integration reports.
- Separate expensive Shiny data construction from app configuration, with cluster
  annotations and metadata presentation settings.

### Fixes

- Align identical gene sets by ID when combining differently ordered inputs for
  join pseudobulks and Shiny data; warn on reordering and reject invalid gene sets.
- Track resolved QC thresholds rather than unrelated metadata/custom-file mtimes,
  avoiding unnecessary threshold regeneration after knee/doublet-option edits.
- Bound Shiny UMAP KDE memory using chunked Gaussian-kernel accumulation while
  preserving the full-data bandwidth, grid and estimate within numerical tolerance.
- Preserve incomplete mapping-knee panel layout and improve workflow dependency,
  environment, metadata and input read-only handling.

### Development and validation

- Release validation on SHPC: 316 Python tests passed, 4 skipped; eight R
  regression scripts passed, covering Shiny alignment/density/label placement,
  scDblFinder options, doublet diagnostics and mapping knee/shin plots. The MkDocs
  build, release version/notes gate and scprocess-shpc skill validation passed.

- Introduce the release policy, rename CLAUDE.md to AGENTS.md and require agents
  to establish version impact and prefer compact commits and clean run checkouts.
- The full 3.35-million-cell Siletti Shiny build passed the earlier KDE memory
  failure but timed out later during sample pseudobulk construction at 45 minutes.
  A completed atlas-scale app build is not claimed by this release.
- The planned Siletti classifier and remaining broad wording cleanup are deferred
  to subsequent commits/releases; no new Siletti model is included here.

Earlier releases: see the immutable GitHub release records and their published tags.
