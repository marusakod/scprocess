# Built-in adult human CNS classifier

`human_cns` is the centrally maintained Siletti supercluster classifier. Normal
`scprocess setup` installs it; no user registration block is needed:

```yaml
label_celltypes:
  - labeller: scprocess
    model: human_cns
```

This applies to project, join and zoom configurations. User-trained models can
still be registered separately in an installation's setup YAML.

## Model identity and installation

The first model is `human_cns-v1`, with 31 classes and 2,134 selected genes.
`resources/xgboost/human_cns.json` pins the download URL, archive checksum and
individual file checksums. Setup verifies the complete archive before installing
its contents and retains the manifest and model card beside the weights.
Predictions and aggregated labels record `model_id` and `model_sha256`.
Replacing weights requires a new immutable model ID and an explicit software
release note; the friendly name `human_cns` is pinned by each scprocess version.

The model download is prepared for a separate `human_cns-v1` release in
`marusakod/scprocessData`. For the 0.3.0 development candidate, publication of
that asset is pending approval from both maintainers. The download URL must be
verified before the scprocess release is approved. Existing installations of
the legacy R-format CNS model need setup rerun for the native JSON model.

## Training and evaluation

The distributed model uses all four Siletti donors, with whole tissue tubes
withheld for early stopping. It samples up to 2,000 training and 500 validation
cells per class, balanced over original cluster, donor and tube combinations.
Original supercluster labels are used without cluster majority voting. Splatter
is included, Misc neuron is excluded, Miscellaneous is split into Leukocyte and
Misc neuron, and Vascular/Fibroblast remain distinct.

Four separate models held out one complete donor each. Test accuracies were
99.44%, 99.17% and 99.40% for donors 2–4, and 96.60% for donor 1, which contains
motor cortex only. These evaluate agreement with publication labels among cells
retained by scprocess QC. The distributed all-donor model has no independent
test donor of its own. Evaluation tables and the model card accompany the weights.

## Input contract and interpretation

Use raw count inputs with the `human_2024` reference identifiers. Inference
normalises each cell using all input genes before selecting model features:
`log1p(10000 * counts / library_size)`. Missing selected genes are filled with
zeros. Reordered input genes are aligned by identifier. Model fitting used
XGBoost 2.0.3 on CUDA; inference runs on CPU in the managed environment.

This is a closed-set classifier of adult human CNS nuclei. Its maximum class
probability is not calibrated confidence, and excluding Misc neuron does not
teach it to reject unknown cell types. Check class support and transfer to
different ages, disease states, preparations and anatomical sampling.
