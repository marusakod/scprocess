# Per-stage iteration variables

## Problem

The pipeline currently uses two variables (`RUN_VAR`, `BATCH_VAR`) to control iteration across all stages. `RUN_VAR` is derived from `demux_type`; `BATCH_VAR` is user-configured. This creates implicit coupling: every stage must reason about `demux_type × batch_var` combinations, leading to branching that is hard to extend. Adding a new data type (e.g. 10x Flex, where FASTQs are per-pool but h5 outputs are per-sample) requires threading new conditionals through the entire pipeline.

## Proposal

Replace `RUN_VAR` / `BATCH_VAR` with per-stage iteration variables:

| Variable | Stage | Meaning |
|---|---|---|
| `MAPPING_VAR` | Alignment | Unit for FASTQ input |
| `AMBIENT_VAR` | Ambient removal + QC | Unit for h5 output / ambient / QC |
| `INTEGRATE_VAR` | Integration onward | Unit for batch correction / clustering |

Note: `QC_VAR` is always equal to `AMBIENT_VAR` (QC operates on whatever the h5 represents), so they share a single variable.

## Values per data type

| Data type | MAPPING | AMBIENT | INTEGRATE |
|---|---|---|---|
| Standard (`demux_type: none`) | sample_id | sample_id | sample_id |
| HTO (sample-level integration) | pool_id | pool_id | sample_id |
| HTO (pool-level integration) | pool_id | pool_id | pool_id |
| Custom demux (sample-level) | pool_id | pool_id | sample_id |
| Custom demux (pool-level) | pool_id | pool_id | pool_id |
| Flex / pre-demuxed (new) | pool_id | sample_id | sample_id |

## Derivation

All three variables are derived from `demux_type` and `int_batch_var` in the config. The user only configures `demux_type` and `int_batch_var`; the per-stage vars are computed:

```python
def get_stage_vars(config):
    demux_type = config['multiplexing']['demux_type']

    if demux_type == "none":
        mapping_var = "sample_id"
        ambient_var = "sample_id"
    elif demux_type in ("hto", "custom"):
        mapping_var = "pool_id"
        ambient_var = "pool_id"
    elif demux_type == "flex":
        mapping_var = "pool_id"
        ambient_var = "sample_id"

    integrate_var = config['integration']['int_batch_var']

    # validate: can't integrate at coarser granularity than QC
    if ambient_var == "sample_id" and integrate_var == "pool_id":
        raise ValueError(
            "cannot integrate at pool level when post-mapping data is already per-sample"
        )

    return mapping_var, ambient_var, integrate_var
```

## Stage handoffs

When adjacent stages have different variables, a mapping is needed. The metadata file provides the pool_id → sample_id relationship.

**Same var across boundary:** output of one stage is directly the input of the next (same wildcard).

**Var becomes finer (pool → sample):** this can happen at two boundaries:
- mapping → ambient (flex case): a "split" rule converts pool-level alignment output into per-sample h5 files
- QC → integration (HTO/custom case): the merge_qc / integration input functions use the pool→sample mapping from metadata

**Var becomes coarser (sample → pool):** this should not happen and is rejected by validation.

## Impact on branching

Each stage only checks its own variable, not `demux_type × batch_var`. For example, `.filter_qc` simplifies from a 6-way branch to:

```r
if (ambient_var == "sample_id") {
    # no demux filtering needed — data is already per-sample
    keep_idx = coldata_in$scdbl_class == 'singlet'
} else if (ambient_var == "pool_id") {
    # demux filtering needed
    keep_idx = (coldata_in$scdbl_class == 'singlet') &
               (coldata_in$demux_class == 'singlet') &
               (coldata_in$valid_hto %in% c(TRUE, NA))
}
```

The demux_type-specific logic stays in `.add_demux_metadata` (which already handles it) and the derivation function above.

## Migration path

1. Add the derivation function. Set `MAPPING_VAR = RUN_VAR`, `AMBIENT_VAR = RUN_VAR`, `INTEGRATE_VAR = BATCH_VAR`. No behavior change.
2. Replace `RUN_VAR` / `BATCH_VAR` references in rules with stage-specific vars. Still no behavior change.
3. Add the flex data type: allow `MAPPING_VAR != AMBIENT_VAR`, add a split rule between mapping and ambient.

## Related: per-stage run lists

Each stage needs its own list of runs/batches:

```python
MAPPING_RUNS   = list(MAPPING_PARAMS.keys())      # pools or samples
AMBIENT_RUNS   = list(AMBIENT_PARAMS.keys())       # samples (flex) or pools (HTO)
INTEGRATE_BATCHES = list(INTEGRATE_PARAMS.keys())   # samples or pools
```

Plus the cross-stage mappings:

```python
MAPPING_TO_AMBIENT   = { pool: [samples] }   # only needed when MAPPING_VAR != AMBIENT_VAR
AMBIENT_TO_INTEGRATE = { pool: [samples] }   # only needed when AMBIENT_VAR != INTEGRATE_VAR
```
