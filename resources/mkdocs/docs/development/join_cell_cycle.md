# GPU Harmony and optional cell-cycle correction in joins

The integration environment pins RAPIDS single-cell 0.16.1 with CUDA 12 and
RAPIDS 25.12. Single or multiple `int_batch_var` covariates run on the GPU when
`int_use_gpu: true`; a theta array follows the covariate order. Covariates remain
independent, not concatenated into a single batch. CPU integration still works
when GPU execution is disabled. GPU Harmony explicitly selects `harmony1` to
retain the original algorithm and the existing five-iteration limit.

Cell-cycle correction is **disabled by default**. To enable it in join or joined-zoom YAML:

```yaml
integration:
  int_pca_method: bpcells
  int_cell_cycle_correction: within_sample
  int_cell_cycle_harmonics: 2
  int_cell_cycle_ridge_lambda: 0.1
```

The last two fields default to 2 and 0.1 only when correction is enabled.
Each retained cell must have a finite source-project `tricycle_theta` score.
Scores are aligned by cell ID from the joint metadata, not recomputed. Missing,
duplicate or degenerate scores stop PCA explicitly. Review phase calibration
across source projects before fitting shared cyclic effects: inconsistent phase
origins are not repaired by join regression. Joins do not need a new `cell_cycle`
configuration block; score cells in the source projects first.

The existing BPCells ridge implementation removes sine/cosine cyclic components
from scaled expression **before PCA**, then Harmony integrates the corrected PCs.
It remains disk-backed and does not materialize a dense cell-by-gene matrix.
Coefficient fitting uses within-sample variation, with project-prefixed sample
IDs. `within_sample` preserves sample means; `global` additionally removes
predicted differences due to phase composition. Correction can remove real
biology, so compare corrected and uncorrected integrations before adopting it.
Counts used for markers, labels and downstream analyses remain unchanged.

PCA exports `cell_cycle_regression_<tag>_<date>.csv.gz` and
`cell_cycle_regression_genes_<tag>_<date>.csv.gz` with model settings, cell-set
fingerprint, coefficients and removed-component variance diagnostics. BPCells
0.3.1 remains required by the validated low-rank residual implementation.
Enabling correction invalidates the PCA and its downstream products; existing
completed results should not be rerun or overwritten without authorization.
Changing the GPU environment causes Snakemake to create a new hashed environment;
it does not modify an old environment in place.

Joined zooms inherit the parent correction settings unless overridden and fit
the model again on their selected cells. They use the same PCA rule and source
scores as ordinary joins. Omission in a zoom inherits its parent setting; it
does not disable correction on a corrected parent. Missing scores fail before
integration, including when only some source projects supplied scores.
