# Changelog

## 0.2.0 — 2026-10-25

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
