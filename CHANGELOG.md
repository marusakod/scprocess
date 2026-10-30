# Changelog

## 0.6.0 — Unreleased development candidate

- Update the shared Python rule environment from Python 3.10 to 3.12,
  including its ABI pin. Setup checksum verification uses hashlib.file_digest;
  Snakemake creates a new environment from the updated specification.

- Replace `cell_cycle.cyc_target_prop` with required `cyc_density_percentile`
  (0–100). The selected cell-density percentile sets the full-weight threshold;
  denser cells are weighted inversely by density. The diagnostic grid is now
  `cyc_density_percentile_grid`, default `[0, 0.1, 0.3, 1, 3, 10, 100]`.
  This changes origin selection: old proportions cannot simply be renamed or
  multiplied by 100. Saved origin/sensitivity tables and reports use percentile
  metadata; regenerate aggregation and affected downstream results after choosing
  a percentile. Per-run projections can be reused. No old-input alias is provided.

- Remove the retired hippocampus classifier download and its separate installer.
- Require canonical `label` input columns and the optional `label_specs_f` /
  `coarse_label_specs_f` tables; remove the old training `label_map_f` interface.
  Coarse evaluation remains separate from model fitting.
- Remove deprecated `user.profile`, `shiny.var_names` and former tricycle/
  join-zoom resource names. Use the current schema names.
- Require saved doublet diagnostics and adaptive/manual integration threshold
  metadata when rendering reports; remove support for older incomplete outputs.
  Recovery guidance for future agent sessions is recorded in AGENTS.md.
- Download support data, reference annotations/indices and tutorial inputs from
  pinned Zenodo records. Verify support/model archives and members before installing.
- Reuse source-project XGBoost calls in joins only when the model checksum
  matches; retain identity fields and reject mixed-weight prediction files.
- Ship the reviewed 32-label/2,386-gene human_cns-v1 replacement, including its
  fine/coarse specifications and evaluation summaries; fitted weights are unchanged.
- Separate support-data updates from Cell Ranger whitelist extraction. Preserve
  other registered custom classifiers when refreshing the installed model list.

## 0.5.0 — Unreleased development candidate

- Add `marker_genes.mkr_dispersion_use_clusters` (default `true`) to projects,
  joins and zooms. Setting it to `false` omits cluster effects only during
  dispersion estimation; filtering, normalization and marker contrasts remain
  unchanged. Fit checkpoints distinguish the two settings.
- Use edgeR.bp as the single marker backend for projects, joins and their zooms.
  Persist BPCells pseudobulks, preparation tables and compact report expression
  data as separate workflow outputs. Existing marker CSV names are retained;
  project/zoom pseudobulk RDS outputs are replaced by `.bpcells` directories.
- Project markers adopt the join calculation: shared TMM offsets, dispersion
  estimation with a cluster design, and streamed one-versus-rest Treat tests.
  Marker statistics may change; this is not a storage-only migration.
- Rank variable pseudobulk genes by TMM-logCPM variance, as joins already do,
  instead of the former project DESeq2 transformation. Reports use this ranking.
- Preserve project doublet exclusions and zoom membership, and aggregate pooled
  libraries by biological sample rather than treating the pool as a sample.
- Replace the obsolete combined marker resource prediction with explicit defaults
  for the split stages; retain fitting memory/time override names.

## 0.4.3 — Unreleased development candidate

- Consolidate doublet-cluster thresholding, streaming diagnostic merging and
  join-report plotting into their existing implementation modules. Keep data
  dependencies and the scientific calculations unchanged.
- Improve processing helper names, function ordering, comments and layout;
  remove unused examples/commented code and place classifier setup helpers
  before CLI dispatch. Setup/download behavior is unchanged.
- Use equivalent base-R normalization calls in PCA, avoiding an undeclared
  magrittr pipe dependency, and distinguish project/join HDF5 helpers.
- Simplify cell-cycle aggregation and its numerical helpers. Move report
  preparation and plotting out of the Rmd into cell_cycle.R; preserve the
  existing report narrative, projections, origin estimates and diagnostic draws.
- Explain classifier evaluation and QC diagnostics for computational biologists.
  Use reference-label terminology, describe training/validation/test groups,
  clarify optional project adjustment and temper absolute ambient-QC claims.
- Validation includes exact threshold and real-cell prediction comparisons,
  unchanged aggregation tables and report PNGs, focused Python/R/workflow
  checks, small report renders and the documentation build. Project marker
  migration and setup/Zenodo migration remain separate future work.

## 0.4.2 — Unreleased development candidate

- Allow `scprocess join ... -r label_celltypes` to run annotation without
  marker genes or the main join report.
- Use installed XGBoost label metadata in joins, mapping both naive and
  majority-vote fine calls to coarse labels when specified. Updating the mapping
  regenerates aggregated labels; fine predictions and voting settings are unchanged.
- Reject an explicitly requested `label_celltypes` step with no classifiers.
- Refactor joined-zoom source preparation and HVG calculation for readability,
  preserving bounded sparse reads, cell selection, outputs and workflow dependencies.
- Remove the unused zoom-name CLI helper.

## 0.4.1 — Unreleased development candidate

- Give every completed version, including patch versions, an immutable annotated
  Git tag. Development-version tags do not publish a release; main promotion and
  release publication still require both maintainers' explicit approval.

- Group fine/coarse diagnostics into nested label/cell-group tabsets and move
  gain/selected-gene plots to the end of the training report.

- Revise XGBoost reports with combined train/validation hex density and local
  validation fractions, square larger UMAPs, log row-proportion confusion
  matrices by cell group, scrollable metrics and legible precision–recall plots
  with full/zoomed axes. Apply the presentation to coarse labels as well.
- Render evaluation tables as HTML and call sampled-out cells "Remaining cells"
  in report prose, tabs and tables. Stored split codes and model results are
  unchanged. Use pretty axis breaks on gain curves.


## 0.4.0 — Unreleased development candidate

### Reviewed training interface

- Consolidate training in `train_xgboost.py`, with Polars tables, explicit
  project/join context, clearer function ordering and documented split units.
- Rename the CLI stages to `prepare`, `train`, `combine` and `specs`; rename
  per-label sampling settings to `*_cells_per_label`, saved plans to
  `*_training_metadata.csv.gz`, and corresponding workflow/resource names.
  Update maintained project, join and zoom callers together.
- Add optional standard/coarse label specification CSVs with independent order
  and named colours. Carry them through reports, prediction and local setup;
  regenerate display/coarse summaries from saved predictions without refitting.
- Place model-path lookup in scprocess utilities and built-in installation in
  setup; simplify registry and model-information helpers.

### Scientific behaviour and validation

- Keep the majority-vote minimum as the named fixed 10-cell method constant.
  Break validation-unit rescue ties by cell count then unit ID, replacing the
  former incidental quicksort tie order; seeded selections can change for ties.
- Require informative training features, required metadata and valid complete
  specification CSV rows; remove silent fallbacks. Gene-type exclusion uses
  native Polars regular expressions.
- Make report gain cutoffs reflect the fitted model's actual target and selected
  gene count. Plot final-model genes only and paginate small models correctly.
- Supply schema defaults for join training preparation and import its label
  specification helper, verified with the complete native join workflow on
  1,440 real Siletti cells from six libraries, including the HTML report.
- Development candidate only. No main promotion, release tags, publication or
  replacement public model asset is included.


## 0.3.0 — Unreleased development candidate

### Training features and interface changes

- Add optional standard/coarse label specification CSVs, independent row ordering
  and named colours, carried through reports and local installation. Standardise
  new label columns as `label`, retaining legacy input names. Derive coarse
  evaluation from saved predictions without refitting or repeating inference.
- Keep optional CLI features omittable and consolidate model registry,
  installation and label-specification helpers into `train_xgboost.py`.
- Resolve current training tags during setup and support direct local source
  directories with an independent installed classifier name.

- Consolidate XGBoost planning, fitting and evaluation into `train_xgboost.py`;
  keep the existing commands, workflow stages and scientific behaviour.
- Use Polars for training plans and evaluation summaries, retaining seeded cell
  selections and model outputs. Pass named settings through project/join/zoom
  rules, with defaults defined in their schemas. Reorder the trainer by stage,
  mark internal helpers with underscores and improve code/rule layout.

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
- Add the built-in human_cns-v1 adult human CNS classifier (31 Siletti labels),
  pinned archive/member checksums, model metadata and donor-held-out evaluations.
- Use the shared XGBoost registry in project, join and zoom validation; remove
  the legacy R-format human_cns lookup. Save model ID/checksum in predictions.
- Development candidate only: model-asset publication and the scprocess release
  remain subject to explicit approval from both maintainers.

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
