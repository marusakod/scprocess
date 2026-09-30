# Joined-zoom implementation

For selection and analysis choices, see [Zooming a completed join](../reference.md#join-zoom).
These notes describe the implementation and its dependency contracts.

## Workflow and identifiers

`get_join_zoom_specs()` validates the subset YAMLs; `get_join_zoom_context()`
resolves inputs and inherited settings. `join_zoom.smk` prepares source tables
and HVGs, then includes `join.smk` for matrix assembly, PCA, integration, marker
testing and reporting. The CLI runs requested zooms sequentially with the
parent's profile and project-level Snakemake locking.

Selection is restricted to cells in the completed parent integration. Exact
cell IDs are preserved, and preparation verifies that each selected cell occurs
in exactly one source H5AD. Temporary per-project tables remove one sample-ID
prefix so the join builders restore the original joined sample IDs. Batch
wildcards use numeric keys rather than parsing IDs on underscores.

Parent clusters and embeddings receive a `parent_` prefix. In particular,
Harmony PCs are named `hmny_pca_*`; their parent counterparts must be
`parent_hmny_pca_*`. Metadata is reattached after integration in the new output
order. Surplus parent PCs must never be appended as new zoom PCs when the zoom
uses fewer dimensions.

## Counts, HVGs and marker models

The source H5ADs are read-only. The subset manifest and link directory include
only batches with selected cells; rebuilding the directory does not change
source files. Empty projects remain valid after selection.

HVG statistics are accumulated within each sample across H5AD files before
fitting the existing Seurat-style VST trend. The supported method is `sample`.
`hvg_chunk_size` defaults to 2000 cells for row-oriented inputs; CSC inputs are
read in bounded gene blocks along their efficient storage axis. Shared PCA
genes must be measured in every retained sample. This chunking applies to HVG
calculation; matrix assembly subsequently uses the shared join builder.

Native BPCells PCA uses the join's library-size normalization, log transform,
centering and scaling. Marker testing uses the join's disk-backed edger.bp path.
If filtering leaves only one project, the constant project effect is omitted
from both preparation and marker fitting.

Native PCA conversion and marker pseudobulking keep HDF5 inputs open read-only
during lazy evaluation. Tests check source modification times and file hashes
as well as numerical output. Readers must not repair an input in place or reset
its timestamps to hide changes.

## Rerun behaviour and resources

Preparation tracks selection, QC, source inputs and the integration settings
it consumes. Marker settings, presentation settings and resource requests do
not themselves invalidate cell selection. `effective_config.yaml` records the
full requested settings through a separate rule; it is not a completion flag.
Actual input changes remain tracked. Older runs may require a one-time rerun
when their dependency records differ.

Rule logs are under `.log/join_zoom/{zoom.name}/` and benchmarks under
`.resources/join_zoom/{zoom.name}/`. The shared `join_pca`, `join_integration`
and `join_marker_genes_calc` resource/profile settings also govern zooms.
Preparation and HVG rules have `gb_join_zoom_*` and `mins_join_zoom_*` overrides.
The Slurm profile controls GPU allocation separately from `int_use_gpu`.

Zoom report helper copies are isolated from the parent's project-level copies.
Shiny app identifiers and cluster annotation paths are not inherited from the
parent. Unsupported subset options are rejected by the joined-zoom schema.
