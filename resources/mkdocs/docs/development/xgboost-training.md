# Training an XGBoost classifier

Use the `train_xgboost` block within a completed project, join or zoom YAML.
Run `scprocess run CONFIG -r train_xgboost` or
`scprocess join CONFIG -r train_xgboost`. Training reuses the existing expression
and integration outputs; review a dry run before launching.

```yaml
train_xgboost:
  xgboost_tag: siletti_superclusters_holdout_donor2
  annots_f: data/configs/cell_labels.csv.gz
  label_specs_f: data/configs/label_specs.csv
  coarse_label_specs_f: data/configs/coarse_label_specs.csv
  excluded_labels: [Misc neuron, Splatter]
  test_group_var: donor_id
  test_groups: [siletti_donor_2]
  validation_group_var: tube_group
  validation_fraction: 0.2
  validation_variable_combinations: [donor_id, project_id]
  train_cells_per_label: 2000
  validation_cells_per_label: 500
  sampling_variable_combinations: [original_cluster, donor_id, tube_group]
  seed: 20261025
```

The annotation CSV needs unique `cell_id` and `label` columns. Grouping and
sampling columns can be supplied there or in the integration table. Shared
columns must agree. `tube_group` should identify a biological tissue tube uniquely
across donors. The workflow inherits processing identifiers from its project or
join context; do not set `batch_var` in the training block.

Metadata preparation, fitting and full-data evaluation are implemented in
`scripts/train_xgboost.py`, through its `prepare`, `train` and
`combine` commands. The `specs` command exports display metadata and
broader-label metrics from saved predictions. They remain separate workflow stages so that
the saved training metadata can be inspected before fitting.


## Label specifications

`label_specs_f` and `coarse_label_specs_f` are optional. All relative paths resolve
from the project folder; `data/configs/` is the suggested location in examples,
not a special search directory. Omitting both settings leaves coarse labels
disabled. Existing files using `annotation` or `fine_label` as their label column
are accepted; new outputs use `label`.

For example, `data/configs/label_specs.csv`:

```csv
label,colour,coarse_label
Astrocyte,#74C476,Astroglia
Bergmann glia,#238B45,Astroglia
Microglia,#A1887F,Immune
Leukocyte,#6D4C41,Immune
```

And `data/configs/coarse_label_specs.csv`:

```csv
coarse_label,coarse_colour
Immune,#A1887F
Astroglia,#74C476
```

Rows in each file independently set its display order. Supply every model label
in the standard table; it can also describe excluded/test-only labels. `colour`
and `coarse_label` may be omitted. Without a coarse mapping, the standard table
controls ordering and colours only. A coarse mapping without a coarse table uses
first-appearance order and the explicit default palette. Colours, when supplied,
must cover the whole corresponding table and use `#RRGGBB` or `#RRGGBBAA`.
Duplicate/missing labels, incomplete mappings and undefined coarse groups fail
validation. The example above illustrates four labels; a real model's table must
include all of its labels.

The model still fits and predicts standard labels. Each coarse label is a direct
lookup from a standard label, applied equally to reference and predicted labels;
there is no voting or probability summation. Both specification tables are saved
beside the model, retained during installation and used by reports. The separate
specification stage reads saved predictions, so changing display metadata does
not require fitting or inference again. Broader-label evaluation CSVs and report
views distinguish independent tests, validation and development holdouts.

Legacy `label_map_f` remains supported, but cannot be combined with
`label_specs_f`. New configurations should use the two specification settings.

## Splitting and sampling

The preparation rule writes `_training_metadata.csv.gz`, `_split_coverage.csv`,
`_label_counts.csv`, and planning metadata before reading expression. The plan
includes the original and training labels, processing batch, split membership,
selection flag and training-class support. Inspect these files before fitting
large references. Changing model hyperparameters reuses the training metadata when its inputs
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
labels below `min_cells_per_label` are excluded from fitting. Test labels absent
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
in the training metadata before enabling it. The Siletti demonstration leaves it disabled.

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

- Rename `ref_tag` to `xgboost_tag`. Prefer `label` over `annotation` in input CSVs.
- Remove `refine_labels`. Omit `cluster_majority_vote_resolution` to disable
  smoothing, or set an explicit saved resolution to enable majority voting.
- Rename `purity_threshold` to `cluster_majority_vote_min_purity` when voting.
- Replace `n_cells_per_type` with separate `train_cells_per_label` and
  `validation_cells_per_label`. Sampling now happens after group assignment.
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

Planning, fitting and specification export pass named settings from the validated
configuration to `train_xgboost.py`. Defaults live in the schemas, including
`validation_group_var: batch` when no grouping variable is supplied. Table
preparation and evaluation use Polars; pandas is limited to AnnData metadata.

The workflow passes `is_join` explicitly to the planning command. Project and
project-zoom training use false; join and joined-zoom training use true. This
controls project metadata and batch-key qualification independently of whether
`batch_var` names samples or pools. It is internal workflow context, not an
additional YAML setting.

The internal CLI requires essential inputs, workflow context and numerical
settings supplied by schema validation. Optional feature arguments can be
omitted: exclusions, test groups, variable combinations, cluster voting and
label mapping/specification files. Missing list arguments disable their feature;
missing file arguments use no mapping/specification. Test grouping still requires
both a variable and a nonempty list of groups. Numerical defaults are not
repeated in argparse.

Model installation, registry identities and label specification helpers also
live in `train_xgboost.py`. Training-only libraries are imported inside the
functions that use them, so setup and configuration checks can reuse these
helpers without loading AnnData or XGBoost.
