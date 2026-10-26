# Training an XGBoost classifier

Use the `train_xgboost` block within a completed project, join or zoom YAML.
Run `scprocess run CONFIG -r train_xgboost` or
`scprocess join CONFIG -r train_xgboost`. Training reuses the existing expression
and integration outputs; review a dry run before launching.

```yaml
train_xgboost:
  xgboost_tag: siletti_superclusters_holdout_donor2
  annots_f: data/classifier_annotations.csv.gz
  excluded_labels: [Misc neuron, Splatter]
  test_group_var: donor_id
  test_groups: [siletti_donor_2]
  validation_group_var: tube_group
  validation_fraction: 0.2
  validation_variable_combinations: [donor_id, project_id]
  train_cells_per_type: 2000
  validation_cells_per_type: 500
  sampling_variable_combinations: [original_cluster, donor_id, tube_group]
  seed: 20261025
```

The annotation CSV needs unique `cell_id` and `annotation` columns. Grouping and
sampling columns can be supplied there or in the integration table. Shared
columns must agree. `tube_group` should identify a biological tissue tube uniquely
across donors. The workflow inherits processing identifiers from its project or
join context; do not set `batch_var` in the training block.

## Splitting and sampling

The planning rule writes `_training_plan.csv.gz`, `_split_coverage.csv`,
`_label_counts.csv`, and planning metadata before reading expression. The plan
includes the original and training labels, processing batch, split membership,
selection flag and training-class support. Inspect these files before fitting
large references. Changing model hyperparameters reuses the plan when its inputs
and planning options are unchanged.

`test_groups` names values of `test_group_var` reserved for independent testing.
Neither their cells nor their labels enter feature selection or early stopping.
Omit both settings for a final candidate using all donors. To evaluate several
held-out donors, use one configuration/output directory per donor; biological
replication comes from donors, not cell counts.

Among the remaining groups, select whole `validation_group_var` units for early
stopping. Without that option the processing batch is used. Within each joint
`validation_variable_combinations` combination, approximately
`validation_fraction` of groups are selected, keeping at least one group in
training. A validation unit cannot cross combinations or the test boundary.
Single-unit combinations stay in training. Rare labels can rescue whole units
from validation to training; rescues and missing validation support are reported.
There is no automatic cell-level fallback.

After splitting, cells are sampled independently within each label and split.
`sampling_variable_combinations` balances joint combinations of its columns.
It is not hierarchical: when combinations outnumber the cap, combinations are
sampled at random. Available cells are redistributed without duplication.
Without sampling variables, sampling is uniform within the label. Test cells
are not downsampled. Small labels retain their available cells; development
labels below `min_cells_per_type` are excluded from fitting. Test labels absent
from training remain evaluation errors rather than being dropped.

Training and validation caps default to 2,000 and 500 cells per label. These are
sampling caps, not a model of natural prevalence or calibrated class priors.
Full-data predictions distinguish selected training/validation cells, independent
test groups, excluded labels and unannotated cells. Unselected development cells
are summarised as development holdouts, separately from independent test results.

## Cluster majority voting

Original labels are preserved when `cluster_majority_vote_resolution` is absent.
Set a saved integration resolution explicitly to enable majority voting:

```yaml
  cluster_majority_vote_resolution: 2
  cluster_majority_vote_min_purity: 0.9
```

A cluster needs at least ten training annotations and the requested majority
fraction. Its training-derived majority label is applied to development labels;
independent test targets remain unchanged. Majority voting can erase genuine
minority cell types within mixed clusters. Inspect the original/training labels
in the plan before enabling it. The Siletti demonstration leaves it disabled.

## Fitting and outputs

Expression loading selects requested H5AD rows in read-only mode and aligns
identical gene sets by identifier when their order differs. Genes must have
nonzero counts in at least `min_cells_expressed` training cells; validation and
test expression do not determine that filter. Library normalisation uses all
input genes before gene filtering: `log1p(10000 * counts / library_size)`.

The existing two-pass XGBoost procedure selects features using first-pass gain,
then fits the selected genes. Each Booster is sliced at its best early-stopping
iteration before feature selection, prediction or export. Outputs retain the
native model JSON, class CSV, selected genes, gain table, selected-cell predictions
and pseudobulks, with a model card recording configuration and checksums.

Full-data evaluation and the HTML report separate early-stopping validation,
independent test performance and development holdouts. Macro F1 and balanced
accuracy average over labels observed in each evaluation set. Predictions of a
training class cannot recover a genuinely unsupported test class, and the latter
still counts as an error. Agreement with supplied annotations is not a guarantee
of biological correctness. Compare per-label support and donor/region coverage.

## Migration to the 0.3 training interface

- Rename `ref_tag` to `xgboost_tag`.
- Remove `refine_labels`. Omit `cluster_majority_vote_resolution` to disable
  smoothing, or set an explicit saved resolution to enable majority voting.
- Rename `purity_threshold` to `cluster_majority_vote_min_purity` when voting.
- Replace `n_cells_per_type` with separate `train_cells_per_type` and
  `validation_cells_per_type`. Sampling now happens after group assignment.
- Grouped tests and joint-combination balancing are optional additions; when
  enabled they intentionally change selected cells and the evaluation question.

The old YAML keys are rejected rather than silently interpreted. Reference
models have independent immutable IDs/checksums. Training does not publish or
register a model, promote main or create a release.

GPU fits (`use_gpu: true`) use a separate managed environment with the explicit
CUDA-enabled XGBoost build. The training job must also receive a Slurm GPU
allocation. A small CUDA probe runs before expression loading and fails if the
build lacks CUDA or execution falls back to CPU. Full-data prediction may run
on CPU; model weights are portable across the two devices. Model cards record
the XGBoost build and verified training device.

## Empty source batches

An excluded or failed batch may have a zero-byte H5AD placeholder. Full-data
prediction writes a header-only CSV for that batch. Aggregation accepts empty
batch outputs and still requires a prediction for every retained integration
cell. Missing files and invalid nonempty H5ADs remain errors.

For the centrally distributed Siletti model, see the [human CNS classifier](../human-cns.md).
