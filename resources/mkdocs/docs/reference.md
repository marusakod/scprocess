# Reference

This guide describes command-line options and YAML settings for setup, individual
projects, joins and zooms. Start with the [tutorials](tutorial.md) for a worked
analysis. Relative paths in a project configuration are resolved from
`project.proj_dir`; join paths are resolved from `join.proj_dir`.
Use `data/configs/` for configuration inputs if convenient; it is not a special
search directory. Omit optional settings unless you need to change them.

## {{scsetup}} { #scprocess-setup }

**Description**: Download support data and the built-in human CNS classifier
from Zenodo, obtain Cell Ranger resources from 10x Genomics, and prepare the
requested reference transcriptomes and probe-set indices for `simpleaf`.
See [Getting started](setup.md) for installation and data-directory setup.

**Parameters**:

* `-n`/`--dry-run` (optional): list the setup steps without running them.
* `-E`/`--extraargs` (optional): quoted Snakemake arguments, for example `-E="--cores 8"`.
* `-c`/`--rangerurl`:  download link for Cell Ranger (v9.0.0 or higher) available on the [10x
  Genomics CellRanger download & installation
  page](https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions);
  required when Cell Ranger resources have not yet been installed.

### Configuration file

Place `scprocess_setup.yaml` in `$SCPROCESS_DATA_DIR`. It specifies the cluster
profile or local CPU allocation, reference genomes and optional custom models
shared across projects. See [data-directory setup](setup.md#scprocess-data-directory-setup)
for the environment variable and initial configuration. For example:

```yaml
user:
  profile_template: slurm_default
  profile_name: slurm_my_cluster
  your_name: Testy McUser
  affiliation: Example institute
  int_use_gpu: false
arvados:
  arv_instance:   instance_name
ref_txomes:
  tenx:
  - name:       human_2024
    decoys:     true
    rrnas:      true
  custom:
  - name:       custom_genome_name
    fasta:      /path/to/genome.fa
    gtf:        /path/to/genes.gtf
    decoys:     true
    mito_str:   "^mt-"
probe_sets:
  tenx:
  - name:       human_v1
xgboost:
  - config:     /path/to/config-training_project.yaml
    name:       my_classifier
```

##### user

* `profile_template`: the name of the bundled HPC profile template to copy from the {{sc}}
  _profiles_ directory into `$SCPROCESS_DATA_DIR/profiles`. Bundled templates include
  `slurm_default`, `lsf_default`, `pbs_default`, `sge_default`, and `htcondor_default`. Exactly one
  of `profile_template` and `local_cores` should be specified.
* `profile_name` (optional): name of the active profile directory under
  `$SCPROCESS_DATA_DIR/profiles`. Defaults to the value of `profile_template`.
* `local_cores`: number of CPU cores available for local execution (see [Snakemake
  documentation](https://snakemake.readthedocs.io/en/v9.8.0/executing/cli.html) for more details).
  Exactly one of `profile_template` and `local_cores` should be specified.
* `your_name` (optional): author's name. If specified it will be used in the configuration file for
  new projects created with the `scprocess newproj -c` command.
* `affiliation` (optional): author's affiliation. If specified it will be used in the configuration
  file for new projects created with the `scprocess newproj -c` command.
* `int_use_gpu` (optional): whether to use GPU acceleration (`RAPIDS-singlecell`) for integration
  and clustering steps. If `false` the value will be used in the configuration file for new projects
  created with the `scprocess newproj -c` command.

##### arvados

* `arv_instance` (optional): the name of the default Arvados instance for the user. If specified it
  will be used in the configuration file for new projects created with the `scprocess newproj -c`
  command.

##### ref_txomes

Prebuilt human and mouse reference transcriptomes from 10x Genomics can be downloaded with
{{scsetup}} by adding `tenx` to the `scprocess_setup.yaml` file. Valid values for names are
`human_2024`, `mouse_2024`, `human_2020`, `mouse_2020`.

Names and specifications for custom references should be listed in the `custom` section of the
`scprocess_setup.yaml` file. For each `custom` genome users have to provide the following
parameters:

* `name`: name to be used for the reference
* `fasta`: path to FASTA file
* `gtf`: path to GTF file
* `mito_str`: regular expression used to identify genes in the mitochondrial genome (example for
  mouse: `"^mt-"`)

Optional parameters for both `tenx` and `custom` references are:

* `decoys`: whether to include decoy sequences when building the mapping index; default `true`.
  Decoys help distinguish transcripts from similar sequences elsewhere in the genome.

Optional parameter for `tenx` references is:

* `rrnas`: whether or not ribosomal RNAs should be included in the reference. If not specified it
  defaults to `true` for all `tenx` genomes.

!!! note "Impact of custom parameters for `tenx` genomes on `scsetup` runtime"

    When configuring `tenx` genomes with their default values, {{scsetup}} will download prebuilt
    indices for `simpleaf`. However, if the default parameters are modified (e.g., setting `rrnas`
    or `decoys` to `false`), {{scsetup}} will build the indices from scratch during execution, which
    will increase the runtime.


!!! info "More about decoys"
    Reads from unannotated genomic regions can resemble transcripts and be
    assigned to the wrong gene. Decoy sequences give the mapper alternatives
    outside the annotated transcriptome, reducing these incorrect assignments.
    Keep `decoys: true` unless you have a reason to omit them. See
    Srivastava et al.[@Srivastava2020-jb] for the mapping approach.

##### probe_sets

Probe set indices for 10x Flex data are prepared with {{scsetup}} by adding a `probe_sets` section
to the `scprocess_setup.yaml` file. Each entry under `tenx` specifies a named probe set to download
and index. Valid values for `name` are `human_v1`, `human_v2`, `mouse_v1`, and `mouse_v2`.


##### xgboost { #setup-xgboost }

Custom `XGBoost` classifiers trained with {{sc}} (using the `train_xgboost` functionality) can be
registered with {{scsetup}} so they are available for cell type annotation in any project. To
register a user-trained classifier, add an `xgboost` section to the `scprocess_setup.yaml` file.
Each entry specifies a classifier to install. There are two ways to reference a trained classifier:

 * use the `config` parameter to point to the project configuration file that was used to train the
   classifier. {{sc}} will automatically locate the training outputs.
 * use the `src_dir` parameter to point directly to the directory containing the trained classifier.
   Also specify `name`, the identifier to use when annotating other projects.

For a local classifier with display metadata, for example:

```yaml
xgboost:
  - name: my_classifier
    src_dir: /path/to/project/output/my_project_train_xgboost
    xgboost_tag: xgboost_my_classifier
```

For installation from `src_dir`, `xgboost_tag` is the full filename prefix of
the trained model; `name` is the identifier used after installation. Alternatively, supply
`config: /path/to/config.yaml` (and an optional `name`) to locate the model from
its training configuration.
Optional `*_label_specs.csv` and `*_coarse_label_specs.csv` are copied with the
weights.

## {{scnew}}

**Description**: Create a new `workflowr` project directory for {{sc}} outputs.

**Parameters**:

* `name` (positional): name of the new `workflowr` project directory.
* `-w`/`--where` (optional): path to the directory where the new project will be created; defaults
  to the current working directory.
* `-s`/`--sub` (optional): if provided, creates `data/fastqs` and `data/metadata` subdirectories
  within the project.
* `-c`/`--config` (optional): generates a template configuration YAML file. If provided, it must be
  followed by either `sc` (single-cell) or `sn` (single-nucleus) to define standard QC thresholds.
  The generated config always includes `tenx_assay_type` in the `project` section. Additional
  modifiers can be appended:
    + `flex`: sets `tenx_assay_type: flex` and includes `probe_set` instead of `ref_txome`. e.g.
      `scprocess newproj project_name -c sc flex`
    + `multiplex`: adds a `multiplexing` section. When combined with `flex`, sets `demux_type:
      flex`; otherwise leaves `demux_type` blank for the user to fill in (e.g. `hto` or `custom`).
      e.g. `scprocess newproj project_name -c sc multiplex` or `scprocess newproj project_name -c sc
      flex multiplex`

## {{scnewjoin}} { #scprocess-newjoin }

**Description**: Create a new `workflowr` project directory for {{scjoin}} outputs.

**Usage**:

```
scprocess newjoin <name> [-w <where>] [-p <config>...]
```

**Parameters**:

* `name` (positional): name of the new `workflowr` project directory.
* `-w`/`--where` (optional): path to the parent directory where the new project will be created;
  defaults to the current working directory.
* `-p`/`--projects` (optional): one or more paths to existing {{sc}} project configuration files.


## {{scknee}} { #scprocess-plotknee }

**Description**: Create an interactive barcode-rank plot. Requires the saved knee-plot data from the
mapping step.

**Parameters**:

* `sample`: sample_id corresponding to the barcode-rank curve.
* `-k`/`--kneefile`: path to the knee plot data file generated by {{sc}}, e.g.
  `output/[short_tag]_mapping/af_[sample_id]/knee_plot_data_[sample_id]_[date_stamp].csv.gz`.
  Exactly one of `--kneefile` and `--configfile` should be specified.
* `-c`/`--configfile`: path to configuration file used for running {{sc}}. Exactly one of
  `--kneefile` and `--configfile` should be specified.


## {{scrun}} { #scprocess-run }

**Description**: Run {{sc}}.

**Parameters**:

* `configfile` (positional): path to a configuration YAML file.
* `-n`/`--dry-run` (optional): list the analysis steps without executing them. Helpful for checking
  input files and parameters.
* `--create-envs` (optional): only create the conda environments needed for the workflow, without
  running any rules.
* `-E`/`--extraargs` (optional): list of additional arguments to pass to `Snakemake`. Refer to
  [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for a
  detailed explanation of available command-line options.
* `-r`/`--rule` (optional): Specifies which rule {{sc}} should run. The options are:
    + `all`: default; includes all [Core pipeline steps](introduction.md#core-pipeline-steps) plus
      `label_celltypes` (if configured).
    + `mapping`: read alignment and quantification.
    + `ambient`: ambient RNA removal (optional) and cell calling.
    + `demux`: sample demultiplexing.
    + `qc`: cell-level QC filtering and doublet identification.
    + `hvg`: calculation of highly variable genes.
    + `cell_cycle`: tricycle projection, phase estimation and the cell-cycle report (requires a
      `cell_cycle` block).
    + `integration`: dimensionality reduction with PCA, optional batch correction with `Harmony`,
      UMAP and clustering.
    + `marker_genes`: marker gene identification and optional gene set enrichment analysis.
    + `label_celltypes`: cell type annotation using a pre-trained classifier.
    + `train_xgboost`: train a custom `XGBoost` classifier from annotated data.
    + `zoom`: subclustering.
    + `shiny`: build an interactive Shiny app from pipeline outputs.
* `--zoom` (optional): only used with `-r shiny`. Build the Shiny app for a zoom (subclustering)
  output. Pass a zoom name (e.g. `--zoom immune_cells`) or `all` to build apps for every zoom
  defined in the config.


### Configuration file

The first example is a minimal poly(A) project using local FASTQs. The second
adds selected QC, annotation, classifier-training and Shiny settings. Replace
the paths and metadata names with your own. Values in the second example are
illustrative choices, rather than a complete list of defaults.

=== "Minimal project"

    ```yaml
    project:
      proj_dir: /path/to/project
      fastq_dir: data/fastqs
      full_tag: example_project
      short_tag: example
      your_name: Testy McUser
      affiliation: Example institute
      date_stamp: "2026-10-30"
      sample_metadata: data/metadata/samples.csv
      tenx_assay_type: poly_a
      ref_txome: human_2024
    ```

=== "With optional analyses"

    ```yaml
    project:
      proj_dir: /path/to/project
      fastq_dir: data/fastqs
      full_tag: example_project
      short_tag: example
      your_name: Testy McUser
      affiliation: Example institute
      date_stamp: "2026-10-30"
      sample_metadata: data/metadata/samples.csv
      tenx_assay_type: poly_a
      ref_txome: human_2024
      metadata_vars: [condition, donor_id]
    qc:
      qc_min_counts: 500
      qc_min_feats: 300
      qc_max_mito: 0.1
    hvg:
      hvg_n_hvgs: 2000
    integration:
      int_pca_method: bpcells
      int_embedding: harmony
      int_batch_var: sample_id
      int_dbl_threshold_method: adaptive
      int_res_ls: [0.1, 0.2, 0.5, 1, 2]
    marker_genes:
      mkr_sel_res: 0.2
      mkr_dispersion_use_clusters: true
    label_celltypes:
      - labeller: scprocess
        model: human_cns
        hi_res_cl: RNA_snn_res.2
    train_xgboost:
      annots_f: data/configs/annotations.csv
      label_specs_f: data/configs/label-specs.csv
      coarse_label_specs_f: data/configs/coarse-label-specs.csv
      xgboost_tag: my_classifier
      validation_group_var: donor_id
      validation_fraction: 0.2
      train_cells_per_label: 2000
      validation_cells_per_label: 500
    shiny:
      app_title: Example project
      email: researcher@example.org
      default_gene: APOE
      n_keep: 30000
      metadata_labels:
        condition: Treatment group
        donor_id: Donor
      metadata_combns:
        - [condition, donor_id]
      cluster_palette: Renoir
      metadata_palettes:
        condition: Set2
    resources:
      retries: 3
      n_run_mapping: 8
    ```

#### Required parameters

##### project

* `proj_dir`: absolute path to `workflowr` project directory created with the {{scnew}} function.
* `fastq_dir`: path to directory containing FASTQ files. Should be absolute or relative to
  `proj_dir`. Exactly one of `fastq_dir` and `arv_uuids` should be specified.
* `arv_uuids`: list of Arvados UUIDs where fastq files are located. Exactly one of `fastq_dir` and
  `arv_uuids` should be specified.
* `arv_instance`: the name of Arvados instance. Required if `arv_uuids` is defined.
* `full_tag`: full project label, used in output file names.
* `short_tag`: abbreviated project label, used in output directory names.
* `your_name`: author’s name, displayed in HTML outputs.
* `affiliation`: author’s affiliation, displayed in HTML outputs.
* `date_stamp`: start date of the analysis, formatted as `"YYYY-MM-DD"`.
* `sample_metadata`: path to CSV file with sample metadata. Should be absolute or relative to
  `proj_dir`. Spaces in column names are not allowed. Only required column is `sample_id`; values in
  `sample_id` should not contain `_R1`/`.R1` and `_R2`/`.R2` strings and should not overlap (a value
  should not be a subset of any other values).
* `ref_txome`: must match one of the values in the `ref_txome` column of `index_parameters.csv`
  (created by {{scsetup}}). Required for polyA data; must not be specified for Flex data (use
  `probe_set` instead).
* `tenx_assay_type`: assay type. Options are `poly_a` and `flex`. When set to `flex`, `probe_set` is
  required and `ref_txome` must not be specified.
* `probe_set`: probe set to use for 10x Flex data. Required for Flex data. Valid values are
  `human_v1`, `human_v2`, `mouse_v1`, and `mouse_v2`. Must not be specified for polyA data (use
  `ref_txome` instead).

#### Optional parameters

##### project

* `tenx_chemistry`: 10x assay configuration. Accepted values are `3v2`, `3v3`, `3v4`, `5v1`, `5v2`,
  `5v3`, `multiome`, `flexv1` and `flexv2`. `multiome` refers only to gene expression data generated
  with the 10x multiome kit (ATACseq data is not supported). For Flex data, chemistry specification
  is redundant as chemistry is derived automatically based on `probe_set`. When not specified for
  polyA data, chemistry is auto-detected from barcode whitelist overlap.
* `metadata_vars`: A list of column names in the `sample_metadata` file to be used for visualizing
  the distribution of cell annotations across identified clusters and regions of the low-dimensional
  embedding.
* `show_arv_uuids`: Whether to display Arvados UUIDs (`arv_uuids`) in the configuration file details
  box on the index page. If `false`, UUIDs are replaced with "not shown". Defaults to `true`.
* `exclude`: List of all samples that should be excluded from the analysis. Samples can be listed
  under `pool_id` (if multiplexed) or `sample_id`.
* `custom_sample_params`: YAML file with optional custom parameters for each pool or sample (custom
  `tenx_chemistry`, custom `mapping`, custom `ambient` and custom `qc` parameters can be specified
  for each sample). For example:

```yaml
pool_id:
  pool_1:
    project:
      tenx_chemistry: 5v2
    mapping:
      knee1: 4000
      shin1: 400
      knee2: 30
      shin2: 5
  pool_2:
    project:
      tenx_chemistry: 5v2
    mapping:
      knee1: 3000
      shin1: 400
      knee2: 30
      shin2: 5
    ambient:
      cb_total_droplets_included: 20000
      cb_learning_rate: 0.001
      cb_posterior_batch_size: 128 # only applicable if cb_version is v0.3.2
sample_id:
  sample_1:
    qc:
      qc_min_counts: 100
```

##### multiplexing

* `demux_type`: `demux_type` options (default is `none`):
    + `none` if experiment is not multiplexed;
    + `hto` if hto-based demultiplexing of samples should be performed with {{sc}};
    + `ocm` if on-chip multiplexing was used. Requires an `ocm_id` column in sample metadata with
      values `OB1`-`OB4`.
    + `custom` if demultiplexing results will be used as input to {{sc}}; or
    + `flex` for 10x Flex data where samples are pooled within a library and demultiplexed using
      probe barcodes. Requires `tenx_assay_type: flex` in the `project` section.
* `fastq_dir`: path to directory containing HTO FASTQ files. Should be absolute or relative to
  `proj_dir`. If `demux_type` is `hto`, exactly one of `fastq_dir` and `arv_uuids` should be
  specified.
* `arv_uuids`: list of Arvados UUIDs where HTO fastq files are located. Expects `arv_instance` to be
  defined. If `demux_type` is `hto`, exactly one of `fastq_dir` and `arv_uuids` should be specified.
* `feature_ref`: path to CSV file with columns `hto_id` and `sequence`. Required if `demux_type` is `hto`.
* `demux_output`: path to CSV file with columns `pool_id`, `sample_id`, `cell_id`. Optional column
  `class` can be added with values `doublet`, `singlet` or `negative`. Required if `demux_type` is
  `custom`.
* `seurat_quantile`: equivalent to the `positive.quantile` argument of the `Seurat::HTODemux`
  function (see [Seurat documentation](https://satijalab.org/seurat/reference/htodemux) for more
  details). Used if `demux_type` is `hto`.

##### ambient

* `ambient_method`: method for ambient RNA removal. Options are:
    + `decontx_background` (default): runs
      [`decontX`](https://bioconductor.org/packages/release/bioc/html/celda.html) with an estimated
      ambient RNA profile as input.
    + `decontx_cluster`: runs `decontX` with a cluster-based background model.
    + `cellbender`: uses [CellBender](https://cellbender.readthedocs.io) for ambient RNA removal.
    + `none`: skips ambient RNA removal;
* `cell_calling`: method for cell calling when `ambient_method` is `none`, `decontx_background`, or
  `decontx_cluster`. Options are `barcodeRanks` (default) and `emptyDrops`.
* `cb_version`: CellBender version to use with `ambient_method: cellbender`;
  default `v0.3.2`. Use this version for new analyses.
* `cb_max_prop_kept`: maximum proportion of droplets, relative to `--total-droplets-included`, that
  `cellbender` can call as cells. Default is `0.9`, meaning samples are excluded if `cellbender`
  calls more than 90% of `--total-droplets-included` droplets as cells. Applicable only if
  `ambient_method` is `cellbender`. For more information about the `--total-droplets-included`
  parameter see [Cellbender
  documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_learning_rate`: Sets the `--learning-rate` `CellBender` parameter to the specified value;
  applicable only if `ambient_method` is `cellbender`. Default value is `0.0001`. For more
  information about this parameter see [Cellbender
  documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_empty_training_fraction`: Sets the `--empty-drop-training-fraction` `CellBender` parameter to
  the specified value; applicable only if `ambient_method` is `cellbender`. Default value is `0.2`.
  Setting this to a lower value (e.g. 0.1 or 0.05) can help if `CellBender` jobs are failing on
  samples with very few cells. For more information about this parameter see [Cellbender
  documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_posterior_batch_size`: Value of the `--posterior-batch-size` parameter; applicable only if
  `ambient_method` is `cellbender` and `cb_version` is `v0.3.2`. For more information about this
  parameter see [Cellbender
  documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html). Default `128`.
* `cb_expected_cells`: forces the `--expected-cells` `Cellbender` parameter to be consistent across
  all samples; applicable only if `ambient_method` is `cellbender`. For more information about this
  parameter see [Cellbender
  documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_total_droplets_included`: forces the `--total-droplets-included` `Cellbender` parameter to be
  consistent across all samples; applicable only if `ambient_method` is `cellbender`. For more
  information about this parameter see [Cellbender
  documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_low_count_threshold`: forces the `--low-count-threshold` `CellBender` parameter to be
  consistent across all samples; applicable only if `ambient_method` is `cellbender`. For more
  information about this parameter see [Cellbender
  documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).


##### qc

* `qc_min_counts`, `qc_max_counts`: minimum and optional maximum UMIs per cell.
  Defaults `500` and `null` (no upper limit).
* `qc_min_feats`, `qc_max_feats`: minimum and optional maximum detected features
  per cell. Defaults `300` and `null`.
* `qc_min_mito`, `qc_max_mito`: lower and upper mitochondrial proportions.
  Defaults `0` and `0.1`. Proportions use values from 0 to 1: `0.1` means 10%.
* `qc_min_splice`, `qc_max_splice`: lower and upper spliced proportions.
  Defaults `0` and `1`, leaving this measure unrestricted. Appropriate choices
  can differ between whole-cell and nucleus data.
* `qc_min_cells`: minimum cells remaining in a sample after QC to retain that
  sample; default `100`.
* `exclude_mito`: exclude mitochondrial genes from the expression matrix
  used for subsequent analysis; default `true`. Mitochondrial content is
  still measured for QC.
* `dbl_clusters`: whether artificial doublets are generated using clusters.
  Default `null` uses the installed scDblFinder default. `false` disables
  cluster-based generation; `true` first clusters the cells and generates
  artificial doublets between clusters. This can help identify mixtures of
  distinct cell types, but small, biologically real populations also need
  careful review.
* `dbl_dbr_sd`: how strongly the expected doublet rate constrains the calling
  threshold. Default `null` uses the installed scDblFinder default. Setting
  `1` removes the expected-rate penalty, leaving the threshold determined by
  classification of artificial doublets. Consider this when the expected
  rate appears unsuitable for a sample; it can substantially change cell
  retention.
* `dbl_min_feats`: minimum detected features for a barcode to enter
  scDblFinder; default `100`.
* `dbl_min_counts`: additional minimum UMIs for entering scDblFinder;
  default `0` disables this filter. Counts exclude rRNA and, when
  `exclude_mito` is true, mitochondrial genes.

  
??? note "Changing the expected-rate setting for selected samples"

    A YAML file supplied through `project.custom_sample_params` can contain:

    ```yaml
    sample_id:
      SAMPLE_01:
        qc:
          dbl_dbr_sd: 1
      SAMPLE_02:
        qc:
          dbl_dbr_sd: null
    ```

##### pb_empties

* `ambient_genes_logfc_thr`: log-fold change (logFC) threshold used to filter the results of the
  edgeR differential expression test comparing empty droplets to cells. Default `0`.
* `ambient_genes_fdr_thr`: false discovery rate (FDR) threshold used to filter the results of the
  edgeR differential expression test comparing empty droplets to cells. Default `0.01`.

##### hvg

* `hvg_method`: how to combine gene variability across cells; default `sample`. Options:
    + `sample` - calculate highly variable genes per sample, then calculate combined ranking across samples;
    + `all` - calculate highly variable genes across all cells in the dataset; and
    + `groups` - calculate highly variable genes for each sample group then calculate combined
      ranking across groups.
* `hvg_metadata_split_var`: if `hvg_method` is `groups`, which variable in `sample_metadata` should
  be used to define sample groups.
* `hvg_n_hvgs`: number of HVGs to use for PCA. Default `2000`.
* `hvg_exclude_ambient_genes`: if `true`, genes enriched in empty droplets relative to cells will be
  excluded from highly variable genes selection. Default `true`.
* `hvg_exclude_from_file`: path to CSV file with genes to be excluded from HVGs. Should be absolute
  or relative to `proj_dir`. File should contain one column, named either `gene_id` or `symbol`.
  Values in the column should all be present in reference genome.
* `hvg_chunk_size`: number of genes to use for each chunked matrix. Default `2000`.


##### integration

* `int_embedding`: use `harmony` (default) to correct for the chosen batch
  variable, or `pca` to cluster without batch correction.
* `int_batch_var`: the variable identifying batches for Harmony; default
  `sample_id`. For HTO/custom multiplexing, `pool_id` is also an option.
  Join analyses can correct for multiple variables, as described in the
  [join configuration](#configuration-file_2).
* `int_theta`: controls how strongly Harmony favours mixing batches;
  default `0.1`. Larger values encourage more mixing.
* `int_n_dims`: number of principal components used to describe expression
  differences; default `50`.
* `int_res_ls`: clustering resolutions to calculate; default
  `[0.1, 0.2, 0.5, 1, 2]`. Larger resolutions generally produce more, smaller
  clusters. Calculating several lets you compare broad cell classes with
  finer subdivision.
* `int_dbl_res`: resolution of this preliminary clustering; default `4`.
  It is separate from the resolutions used to examine the final cell populations.
* `int_dbl_threshold_method`: `adaptive` (default) estimates the cutoff from
  the observed cluster doublet proportions. `manual` uses your chosen cutoff.
* `int_dbl_cl_prop`: cutoff for the manual method, expressed as a fraction.
  For example, `0.3` excludes clusters with **more than 30%** called doublets.
  It is required with `manual`, has no default, and must be omitted with
  `adaptive`.
* `int_umap_n_iter`: number of steps used to optimize the UMAP layout;
  default `2000`. Applies to project, join and zoom analyses on CPU or GPU.
  Projects use it for both preliminary and final UMAPs. Changing it requires
  recalculating the embedding and affected downstream outputs, without
  changing the requested PCA, Harmony or clustering settings.
* `int_use_paga`: use relationships between clusters to initialize the UMAP
  through partition-based graph abstraction (PAGA); default `true`.
* `int_paga_cl_res`: clustering resolution used for PAGA; default `2` for
  projects, joins and zooms. It must be present in `int_res_ls`. Omit this
  setting if `int_use_paga` is `false`. It is independent of the marker
  resolution, which defaults to `0.2`.
* `int_pca_method`: `bpcells` (default) reads normalized and scaled
  expression from disk, avoiding the need to hold the full matrix in memory.
  `scanpy` calculates PCA in memory on CPU or GPU.
* `int_use_gpu`: use GPU acceleration for supported integration and
  clustering calculations; default `true`. A usable GPU allocation and
  software environment are required. If a GPU is requested but unavailable,
  execution stops with an error; use `false` to request CPU processing.
  Joins with multiple Harmony batch variables use CPU Harmony with
  GPU-accelerated neighbours, clustering and UMAP, and report this choice.
* `int_sce_outs`: also save H5AD outputs as RDS files containing
  `SingleCellExperiment` objects; default `false`.
* `int_cell_cycle_correction` has two options: `within_sample`: removes cell-cycle-associated differences between cells
  within each sample, preserving each sample's average expression on the scale
  used for PCA. Start here when you want to reduce cell-cycle structure before
  Harmony integration. Differences between samples with more or fewer cycling
  cells are retained at this step. `global`: also removes the predicted differences between sample averages
  attributable to their different cell-cycle composition. Consider this when
  those differences dominate the leading PCs. It does not remove all
  differences between samples.
* `int_cell_cycle_harmonics`: how flexibly the correction follows changes
  around the cycle; `1` or `2`, default `2`.
* `int_cell_cycle_ridge_lambda`: limits the strength of the fitted
  correction; default `0.1`. Larger values reduce the correction.


  
##### cell_cycle

* `cyc_density_percentile` sets the density threshold for full weight and must
  be supplied, from `0` to `100`. For example, `1.3` uses the 1.3rd percentile
  of density across cells. Cells at or below that density receive full weight;
  cells in denser areas receive less weight. Tied densities can make the
  percentage at full weight larger than requested.
  Smaller percentiles give less populated areas more influence. At `0`, weights
  are proportional to inverse density; at `100`, every cell has equal weight.
  This is **not** the percentage of cycling cells or a filter that removes cells.
* `cyc_origin_manual` optionally sets the centre yourself, as `[PC1, PC2]` using
  the coordinates in the report. Keep `cyc_density_percentile` so that the estimated
  alternatives are still shown for comparison.
* `cyc_density_percentile_grid`: alternative density percentiles shown in
  the report; default `[0, 0.1, 0.3, 1, 3, 10, 100]`. Your selected value is
  always included. Values are already on a percentage scale; use `1.3` for
  the 1.3rd percentile, rather than `0.013`.
* `cyc_bandwidth_multiplier`: how broadly neighbouring cells are grouped
  when estimating the density of cells in the plot; default `1`. Larger
  values smooth over a wider area and can change the estimated centre.
* `cyc_kde_grid_size`: resolution of the grid used to estimate cell density;
  default `500` bins per axis.
* `cyc_seed`: random seed for the subset displayed in the diagnostic plot;
  default `20230308`. It does not change the centre or cell-cycle positions.


##### marker_genes

* `mkr_sel_res`: clustering resolution used for marker analysis; default `0.2`.
  Choose one of the resolutions already calculated in `integration.int_res_ls`.
* `mkr_min_cl_size`: minimum total cells in a cluster to calculate its markers;
  default `100`.
* `mkr_min_cells`: minimum cells in each cluster–sample pseudobulk for it to
  contribute to marker testing; default `10`.
* `mkr_min_cpm_mkr`: minimum expression for selecting genes shown in marker
  plots; default `50` counts per million (CPM). This is an expression filter
  for the report, not a significance cutoff for the marker test.
* `mkr_not_ok_re`: pattern matching gene types to omit from marker-expression
  plots; default `"(lincRNA|lncRNA|pseudogene|antisense)"`.
* `mkr_do_gsea`: run GSEA; default `true`.
* `mkr_gsea_cut`: FDR cutoff used for GSEA results; default `0.1`.
* `mkr_min_cpm_go`: minimum expression for genes included in GSEA;
  default `1` CPM.
* `mkr_max_zero_p`: maximum proportion of a cluster's pseudobulks with
  zero counts when filtering genes for enrichment; default `0.5`. The zero
  filter is not applied to genes with negative `logFC`.
* `mkr_gsea_var`: how genes are ranked. The default, `z_score`, uses a
  z-score based on signed `log10(FDR)`, combining the direction and
  statistical support of the difference. `logFC` ranks by fold change.
* `mkr_custom_genesets` supplies named gene panels for the report. Each entry
    contains:
    * `name`: a unique name for the panel.
    * `file`: path to a CSV containing `label` (the marker category), and
      `symbol` and/or `ensembl_id`. If the path is omitted, {{sc}} looks for
      `$SCPROCESS_DATA_DIR/marker_genes/{name}.csv`.
* `mkr_dispersion_use_clusters`: include cluster membership when estimating
  variation between samples; default `true`. This separates systematic
  expression differences between clusters from the variation being
  estimated. `false` omits cluster membership from dispersion estimation.
  Gene filtering, TMM normalization and the cluster-versus-rest tests keep
  their existing designs. Changing this setting requires refitting markers
  and updating reports, but not rebuilding counts or normalization tables.
* `mkr_adjust_project_id`: joins only; account for project-level differences
  in both dispersion estimation and marker comparisons. Default `false`.
  `true` adds `project_id` to these models and takes longer to fit. This
  choice is independent of `mkr_dispersion_use_clusters`.

##### label_celltypes

The report shows both per-cell predictions (naive labels) and labels assigned
by the most frequent prediction in each high-resolution cluster (aggregated
labels). Aggregation can reduce noisy individual calls, but can also obscure
rare populations within a cluster; inspect both views.

* `labeller`: specifies the method to annotate cell types; options include:
    + `celltypist` uses one of the models available in
      [`CellTypist`](https://www.celltypist.org/models) for annotation.
    + `scprocess`: use an `XGBoost` classifier for cell type annotation.
* `model`: determines the model to be used based on the selected `labeller`. For list of all
  available `CellTypist` models see `$SCPROCESS_DATA_DIR/celltypist/celltypist_models.csv`. If
  `labeller` is set to `scprocess` the value should match the name of a classifier registered with
  {{scsetup}} (see `$SCPROCESS_DATA_DIR/xgboost/xgboost_models.csv` for a list of available
  classifiers).
* `hi_res_cl`: name of a column containing high-resolution clustering results. It must follow the
  pattern `"RNA_snn_res.n"` where `n` should be replaced with one of the values in `int_res_ls`.
  Default is `"RNA_snn_res.2"`.
* `min_cl_prop`: minimum fraction of a cluster sharing its most frequent predicted label for that
  label to be assigned to the whole cluster. Default `0.5`; clusters below this threshold are
  labelled ambiguous.
* `min_cl_size`: minimum number of cells in a cluster required for that cluster to be labeled. Default `100`.
* `label_specs_f`: one row per label, in the desired display order. Required
  column: `label`. Optional columns: `colour` and `coarse_label`. Should be
  absolute or relative to `proj_dir`. If not specified, {{sc}} will use the
  default label specification included with the classifier (if available).
* `coarse_label_specs_f`: one row per coarse group, in its display order.
  Required column: `coarse_label`; optional column: `coarse_colour`. This file
  requires a coarse mapping in `label_specs_f`.
* `save_cluster_names_file`: if `true`, generate a `cluster_names_for_shiny_*.csv` file mapping
  clusters (at the `marker_genes:mkr_sel_res` resolution) to predicted cell type names. This file
  can be used as the `annotation_csv` in the shiny config. Requires a `marker_genes` block to be
  configured. Default is `false`.


##### train_xgboost

!!! info "Splitting input dataset into training, validation, and test sets"
    
    The partitioning of cells into training, validation, and test sets occurs at the group level i.e. all cells belonging to a specific group are assigned together to a single set.
    
    By default, groups correspond to individual samples (`sample_id`) or pools (`pool_id`), but they can also be mapped to any other metadata variable (such as a batch variable). For example, a metadata variable can optionally be used as `test_group_var`. In this case, all samples corresponding to `test_groups` serve as the test set, held out completely from model training and reserved for diagnostic plots and evaluation.
    
    Similarly, a metadata variable can be specified as `validation_group_var`. A designated fraction of its unique values (`validation_fraction`) is randomly selected, and all samples corresponding to those values form the validation set. Unlike test cells, validation cells directly influence model training (e.g., determining early stopping during boosting).
    
    Any remaining cells not allocated to the test or validation sets are used for training, with up to `train_cells_per_label` cells sampled per class.

* `annots_f` (required): CSV file with a unique `cell_id` column and a `label`
  column containing training labels. Can also supply other metadata columns (any
  values shared with sample metadata must agree). Should be absolute or relative
  to `proj_dir`.
* `excluded_labels`: labels to leave out of training.
* `min_cells_per_label`: minimum number of cells needed
  to retain a label for training; default `20`. Labels with fewer cells are excluded.
* `test_group_var`: metadata column corresponding to the variable with values that should be completely held out during training and validaton; should be supplied together with `test_groups`.
* `test_groups`: values of `test_group_var`.
* `validation_group_var`: metadata column defining the groups to keep together (validation data split is gonna be defined on these groups)
  for validation. By default it uses the processing batch (batch_var?) inherited from the project or join.
* `validation_fraction`: fraction of the groups assigned to validation set; default `0.2`. With ten samples, this initially selects two
  samples (The final cell fraction depends on how many cells each sample contributes)
* `validation_variable_combinations`: metadata columns (other than `validation_group_var` itself) whose combinations should each be represented
  on both sides of the split. Without this, `validation_fraction` draws validation groups from the whole pool at once, which can by chance leave a
  whole condition — e.g. a region or condition that only a few groups elong to — entirely out of validation. Listing `[region, condition]` instead
  draws `validation_fraction` of the groups separately within each distinct `(region, condition)` combination, so every combination contributes to both
  train and validation. Omit for a single, unstratified draw over all groups.
* `sampling_variable_combinations`: balance cell selection within each label across combinations such as `[donor_id, region]`. This happens after group
  assignment, separately for training and validation. Omit to give cells within each label and assigned split an equal chance of selection.
* `train_cells_per_label`: maximum cells to select per label after assigning groups. default: 2000
* `validation_cells_per_label`: maximum cells to select per label after assigning groups; default: 500.
* `cluster_majority_vote_resolution`: an already calculated clustering
  resolution at which to take each cluster's most frequent training label. Omit to keep the original labels.
* `cluster_majority_vote_min_purity`: minimum fraction of labelled training cells that must agree before their label can be assigned to the cluster;
  default `0.9`. For example, at least 90% must share the most frequent label. At least ten training annotations are also needed per cluster — a fixed
  minimum, not a configurable setting. Votes are determined from training cells only, keeping validation and test labels out of that decision.
* `min_cells_expressed`: a gene must be detected in at least this many training cells; default `10`.
* `gene_exclude_re`: pattern matching gene types to exclude; default `"(lincRNA|lncRNA|pseudogene|antisense)"`.
* `gain_threshold`: select genes accounting for this fraction of total gain in the first fit; default `0.9` (90%).
* `min_genes`: lower limit on the selected gene count; default: 100
* `max_genes`: upper limit on the selected gene count; default: 3000
* `seed`: random seed for reproducibility; default `42`.
* `use_gpu`: use a GPU for training; default `false`. 
* `xgboost_tag`: name used in model filenames. If omitted, the file prefix is `xgboost_{full_tag}` for a
  project or `xgboost_{join.name}` for a join.
* `pass1_subsample`: fraction of training cells sampled for each round in the first fit. Default `0.632`. Reused for the second fit as well — there is no
  separate `pass2_subsample`.
* `pass1_colsample_bytree`: fraction of genes considered for each tree in the
  first fit. Default `0.1`.
* `pass1_learning_rate`: how much each added tree changes predictions in the
  first fit. Default `0.1`.
* `pass1_nrounds`: maximum rounds in the first fit. Default `300`.
* `pass1_early_stopping`: stop the first fit after this many rounds without
  improved validation performance. Default `10`.
* `pass2_colsample_bytree`: fraction of genes considered for each tree in the
  final fit, over just the genes selected after the first fit. Default `0.5`.
* `pass2_learning_rate`: how much each added tree changes predictions in the
  final fit. Default `0.05`.
* `pass2_nrounds`: maximum rounds in the final fit. Default `500`.
* `pass2_early_stopping`: stop the final fit after this many rounds without
  improved validation performance. Default `10`.
* `label_specs_f` (optional): one row per label, in the desired display order. Required column: `label`. Optional columns: `colour` and `coarse_label`.
* `coarse_label_specs_f`: one row per coarse group, in its display order. Required column: `coarse_label`; optional column: `coarse_colour`. Requires a
  coarse mapping in `label_specs_f`.


##### shiny

* `app_tag`: identifier used in the app directory, data filenames and build log. It may contain
  letters, numbers, `.`, `_`, and `-`. Defaults to `project.short_tag` for a standard project and
  `join.name` for a join. Main apps are written to `public/shiny_{app_tag}`.
* `app_title`: title displayed in the Shiny app header. Defaults to `short_tag`.
* `email`: contact email shown in the app footer.
* `keyword`: short word used in plot axis labels and descriptions (e.g. `"cells"`, `"nuclei"`).
  Default is `"cells"`.
* `default_gene`: gene symbol displayed by default in the "Explore Genes" tab.
* `n_keep`: number of cells retained in the subsampled UMAP shown in the app. Default is 30000.
* `metadata_vars`: metadata variables to show in the app. If not specified, uses the project-level
  `metadata_vars`. Use this to show a different set of variables in the Shiny app than in the main
  pipeline.
* `metadata_labels`: display labels for metadata variables. Keys are column names, values are
  display labels. Variables not listed keep their column name as the label. Example:
    ```yaml
    shiny:
      app_tag: test_project
      metadata_labels:
        brainregion: "brain region"
        condition: "treatment group"
    ```
* `metadata_combns`: list of metadata variable pairs to display as combined groupings. Each element
  should be a two-element list of variable names from `metadata_vars` (use column names, not display
  labels).
* `home_md`: path to a Markdown file used as the landing page content. Should be absolute or
  relative to `proj_dir`.
* `annotation_csv`: path to a CSV file with columns `cluster`, `cluster_name`, and optionally
  `colour`, defining display names, order, and colours for clusters. Absolute or relative to
  `proj_dir`. During app configuration, its cluster values must exactly match the full set of
  clusters in the built Shiny data, including rare clusters without marker-gene results. Missing or
  extra clusters, duplicated cluster IDs or display names, and missing or blank cluster IDs or names
  cause configuration to stop with an error; this prevents annotations from a different clustering
  resolution or analysis run being applied silently. UMAP subsampling retains at least one cell from
  every cluster.
* `cluster_palette`: name of colour palette applied to clusters when `annotation_csv` is not
  provided. Accepts any name from the Supported colour palettes list.
* `metadata_palettes`: per-variable colour palette configuration. Each key is a metadata variable
  name; the value can be:
    - a palette name (string)
    - an explicit colour list (array)
    - an object with optional `palette`, `colours`, and `values` keys (the `values` list controls
      level ordering; defaults to frequency order)

??? info "Supported colour palettes"

    The following palette names are accepted anywhere a `cluster_palette` or `metadata_palettes`
    entry expects a name:

    | Source | Names |
    |--------|-------|
    | **scprocess built-in** | `nice_cols` (42 colours) |
    | **MetBrewer** | `Archambault`, `Austria`, `Benedictus`, `Cassatt1`, `Cassatt2`, `Cross`, `Degas`, `Demuth`, `Derain`, `Egypt`, `Gauguin`, `Greek`, `Hiroshige`, `Hokusai1`, `Hokusai2`, `Hokusai3`, `Homer1`, `Homer2`, `Ingres`, `Isfahan1`, `Isfahan2`, `Java`, `Johnson`, `Juarez`, `Kandinsky`, `Klimt`, `Lakota`, `Manet`, `Monet`, `Moreau`, `Morgenstern`, `Nattier`, `Navajo`, `NewKingdom`, `Nizami`, `OKeeffe1`, `OKeeffe2`, `Paquin`, `Peru1`, `Peru2`, `Pillement`, `Pissaro`, `Redon`, `Renoir`, `Signac`, `Tam`, `Tara`, `Thomas`, `Tiepolo`, `Troy`, `Tsimshian`, `VanGogh1`, `VanGogh2`, `VanGogh3`, `Veronese`, `Wissing` |
    | **RColorBrewer** | Qualitative: `Accent`, `Dark2`, `Paired`, `Pastel1`, `Pastel2`, `Set1`, `Set2`, `Set3`. Sequential/diverging: `Blues`, `BuGn`, `BuPu`, `GnBu`, `Greens`, `Greys`, `Oranges`, `OrRd`, `PuBu`, `PuBuGn`, `PuRd`, `Purples`, `RdPu`, `Reds`, `YlGn`, `YlGnBu`, `YlOrBr`, `YlOrRd`, `BrBG`, `PiYG`, `PRGn`, `PuOr`, `RdBu`, `RdGy`, `RdYlBu`, `RdYlGn`, `Spectral`. |
    | **ggsci** | `npg`, `aaas`, `nejm`, `lancet`, `jama`, `jco`, `ucscgb`, `d3` / `d3_20`, `d3_10`, `d3_20b`, `d3_20c`, `igv`, `cosmic`, `simpsons`, `rickandmorty`, `futurama`, `tron`, `startrek`, `uchicago`, `frontiers`, `flatui`, `bootstrap` |


##### zoom

List zoom YAML files in the project's `zoom` section to reanalyse selected cells.
HVG, integration, marker, Shiny and XGBoost settings inherit from the parent
project unless overridden. Zoom-specific cell-level QC thresholds are omitted
by default; specifying them adds filtering to the already retained cells.
See [Zooming a completed join](#join-zoom) for the corresponding join analysis,
which has different supported settings and outputs.

Additional parameters include:

* `name`: name of cell subset to be analysed.
* `labels_source`: specifies how a cell subset is defined (required). Options include:
    - `scprocess`: labels assigned by the `XGBoost` classifier (using rule `label_celltypes`)
    - `celltypist`: labels assigned by `CellTypist`(using the rule `label_celltypes`)
    - `clusters`: labels based on clustering results obtained with {{sc}}
    - `custom`: user-defined cell type annotations
* `model`: required if `labels_source` is `scprocess` or `celltypist`.
* `sel_labels`: a list of all labels that define cell types/clusters to be included in subclustering
  (required).
* `labels_col`: name of the column containing the labels used to select cells (required).
* `save_subset_sces`: whether to create `SingleCellExperiment` objects containing cells that have
  been assigned one of the values in `sel_labels`; default is `false`.
* `save_subset_anndata`: whether to create H5AD files containing cells that have been assigned one
  of the values in `sel_labels`; default is `true`.
* `save_cluster_names_file`: whether to generate a resolution-specific
  `cluster_names_for_shiny_*.csv` for the zoom clusters. Names are assigned by majority vote from
  the retained cells' `predicted_label_naive` values and made unique when the same predicted subtype
  is assigned to multiple clusters. This requires a `celltypist` or `scprocess` labels source and
  can be supplied as `shiny.annotation_csv`. Default is `false`.
* `custom_labels_f`: required if `labels_source` is set to `custom`; path to CSV file with columns
  `sample_id`, `cell_id` and `label`.
* `exclude`: optionally excludes batches from this zoom only. Use `sample_id` or `pool_id` to match
  the project's batch variable, for example `exclude: {sample_id: [sample_a, sample_b]}`. This does
  not exclude the batches from the standard analysis or other zooms.

If a value in `sel_labels` is absent from the selected labels column, {{sc}} reports a warning and
continues with the labels that are present. At least one selected label must be present, and the
retained subset must contain enough cells for downstream QC, HVG detection, and integration.

Zoom configs also support optional cell-level QC thresholds to apply stricter filtering to cells
that already passed the main pipeline's QC. When set, cells in the zoom subset that do not meet
these thresholds are removed before HVG detection and integration. All thresholds default to `null`
(no additional filtering):

The zoom report always shows QC distributions. Without zoom-specific thresholds it contains one
distributions view; when thresholds are configured it also shows the before/after filtering views.

When a zoom is defined from `celltypist` or `scprocess` labels, the zoom report also reuses the
retained cells' naive `label_celltypes` predictions. It re-aggregates those predictions by majority
vote against the zoom analysis high-resolution clusters, using the matching parent `label_celltypes`
entry's `hi_res_cl` (default `RNA_snn_res.2`) and `min_cl_prop`. The report shows naive predicted
labels against those high-resolution zoom clusters, zoom-aggregated labels over the zoom UMAP, and
naive-versus-aggregated label totals. Coarse-label panels are included when the source label output
contains coarse predictions. The zoom workflow does not rerun the classifier.

* `qc_min_counts`: minimum UMI counts per cell. Only effective if stricter than the main pipeline
  threshold. Default: not set.
* `qc_max_counts`: maximum UMI counts per cell. Default: not set.
* `qc_min_feats`: minimum detected genes per cell. Default: not set.
* `qc_max_feats`: maximum detected genes per cell. Default: not set.
* `qc_max_mito`: maximum mitochondrial proportion (0-1). Default: not set.
* `qc_min_mito`: minimum mitochondrial proportion (0-1). Default: not set.
* `qc_max_splice`: maximum spliced proportion (0-1). Default: not set.
* `qc_min_splice`: minimum spliced proportion (0-1). Default: not set.

```yaml
zoom:
  name: oligos_opcs
  labels_source: celltypist
  model: Mouse_Whole_Brain
  sel_labels: ["327 Oligo NN", "326 OPC NN"]
  labels_col: predicted_label_agg
  save_subset_sces: true
  save_subset_anndata: true
  save_cluster_names_file: true
qc:
  qc_min_cells: 100
  qc_max_mito: 0.05
  qc_min_counts: 1000
hvg:
  hvg_method: all
shiny:
  app_title: Oligos & OPCs
  n_keep: 20000
  default_gene: Mog
```

##### resources

Use this section to change CPU, memory or runtime requests for individual
workflow steps, for example after a job runs out of memory or reaches its time
limit. Resource settings request an allocation; extra CPUs help only when the
calculation can use them.

* `gb_[rule_name]`: memory requested for the entire job in GB, not per thread.
* `mins_[rule_name]`: maximum runtime requested in minutes.
* `n_[rule_name]`: CPU threads, for rules that expose this setting.
* `retries`: number of retries after failure; default `3`. Requested memory is
  multiplied by `1.5 ** (attempt - 1)`, where the first attempt is `1`.
* `n_run_mapping`, `n_run_mapping_flex`: mapping threads for poly(A) and Flex
  data, respectively; both default to `8`.

For example, to request 16 GB and one hour for project marker fitting:

```yaml
resources:
  n_run_marker_genes: 8
  gb_run_marker_genes: 16
  mins_run_marker_genes: 60
```

By default, many rules estimate resources from benchmarked input sizes or batch
counts. The estimate is `max(floor, intercept + slope × input_size)`, with a
rule-specific input measure. Five per-sample rules have separate estimates for
10x v3 and v4 chemistry. Explicit settings override these estimates; rules
without a fitted estimate use their schema defaults.

??? note "Marker-gene resource names"

    The memory and runtime prefixes above apply to the following stages:

    | Stage | Project rule | Project-zoom rule | Join rule |
    |-------|--------------|-------------------|-----------|
    | Build pseudobulks | `make_marker_pseudobulks` | `zoom_make_marker_pseudobulks` | `join_make_pseudobulks` |
    | Filter and normalize | `prepare_marker_pseudobulks` | `zoom_prepare_marker_pseudobulks` | `join_marker_genes_prep` |
    | Fit marker tests | `run_marker_genes` | `zoom_run_marker_genes` | `join_marker_genes_calc` |
    | Rank variable genes | `calc_marker_hvgs` | `zoom_calc_marker_hvgs` | `join_calc_hvgs` |
    | Prepare expression plots | `prepare_marker_plot_data` | `zoom_prepare_marker_plot_data` | `join_marker_plot_data` |

    Fitting CPUs are controlled by `n_run_marker_genes` for projects and project
    zooms, or `n_join_marker_genes_calc` for joins. Both default to `8`.
    Joined zooms use the join marker settings.

For the complete set of supported resource keys, consult
`resources/schemas/config.schema.json`, `join.schema.json` or `zoom.schema.json`
in your installed {{sc}} checkout. Use current rule names in cluster-profile
overrides too. A dry run shows which rules are needed; inspect child-job logs
to confirm the allocations actually requested.

## {{scjoin}} { #scprocess-join }

**Description**: Integrate multiple completed {{sc}} projects.

**Parameters**:

* `configfile` (positional): path to a configuration YAML file.
* `-r`/`--rule` (optional): `all` (default), `label_celltypes` (annotate cells), `shiny` (build the
  Shiny app from join outputs), `train_xgboost` (train an XGBoost classifier from join outputs), or
  `zoom` (reanalyse a selected population from a completed join).
* `--zoom` (optional): with `-r zoom` or `-r shiny`, select a named joined zoom, or `all`. Without
  this option, `-r zoom` runs every listed zoom; `-r shiny` builds the parent join's app.
* `-n`/`--dry-run` (optional): list the join steps without executing them. Helpful for checking
  input files and parameters.
* `-E`/`--extraargs` (optional): list of additional arguments to pass to `Snakemake`. Refer to
  [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli).
* `--unlock` (optional): unlock the directory if a previous run was interrupted.

### Zooming a completed join { #join-zoom }

A joined zoom reanalyses a selected population from a completed join to reveal
subtypes that may be obscured in the full atlas. It recalculates HVGs, PCA,
Harmony, clusters, UMAP and marker genes. Ambient correction and doublet detection
are not repeated; cells excluded from the parent analysis cannot be restored.

Add a list of zoom YAML files to the parent join configuration:

```yaml
zoom:
  - data/configs/zoom-neurons.yaml
```

For example, `data/configs/zoom-neurons.yaml` could contain:

```yaml
zoom:
  name: neurons
  labels_source: clusters
  labels_col: RNA_snn_res.0.5
  sel_labels: [cl1, cl3]  # replace with labels present in your joined atlas
qc:
  qc_min_cells: 10
hvg:
  hvg_n_hvgs: 3000
integration:
  int_n_dims: 50
  int_res_ls: [0.2, 0.5, 1]
  int_paga_cl_res: 1
  int_umap_n_iter: 2000
marker_genes:
  mkr_sel_res: 0.5
```

Check the dry run, then run the analysis from your usual compute environment:

```bash
scprocess join data/configs/join-atlas.yaml -r zoom --zoom neurons -n
scprocess join data/configs/join-atlas.yaml -r zoom --zoom neurons
scprocess join data/configs/join-atlas.yaml -r shiny --zoom neurons -n
scprocess join data/configs/join-atlas.yaml -r shiny --zoom neurons
```

Omit `--zoom neurons` to run all listed zooms sequentially. Running the parent
join does not run its zooms. Each zoom has separate outputs and a separate Shiny
app; adding one does not itself require rerunning the parent analysis.

**Selecting cells**

* Supply the four selection fields shown above. Use a unique `zoom.name` starting
  with a letter or digit, followed by letters, digits, `.`, `_` or `-`.
* Select parent clusters with `labels_source: clusters`, or existing `celltypist`
  or `scprocess` labels with the corresponding `model`. For custom annotations,
  use `labels_source: custom` and `custom_labels_f`: a CSV/CSV.gz containing unique
  `cell_id` values and the chosen labels column. A file of cell IDs alone is not
  sufficient. Partial annotations are allowed, but requested labels must exist.
* `zoom.exclude.sample_id` excludes exact joined sample IDs, including their
  project prefixes. A zoom may retain cells from only one source project.
* Optional `qc` thresholds further filter cells using the parent QC metrics.
  `qc_min_cells` defaults to 10 selected cells per sample; at least two are always
  required for variance estimation. Check `selection_summary.csv` for loss of
  samples that could affect biological comparisons.

**Choosing analysis parameters**

Settings inherit from the parent unless overridden. In particular:

* HVGs are ranked within samples using the selected cells' clean counts
  (`hvg_method: sample`). `hvg_n_hvgs` sets how many enter PCA. To exclude genes,
  set `hvg_exclude_from_file` to a CSV with gene IDs in its first column; ambient
  gene exclusions are not automatically inherited.
* Choose `int_batch_var` for the technical variation Harmony should remove.
  It must be populated for every selected cell in the parent's cell or sample
  metadata. Avoid adjusting away the biological distinction you want to study.
* Consider fewer PCs (`int_n_dims`) for small populations. UMAP uses 2000
  iterations unless changed in the parent or zoom. Choose `mkr_sel_res`, and
  `int_paga_cl_res` when PAGA is enabled, from the calculated `int_res_ls` values.
* Markers use raw-count pseudobulks from the new clusters. Project adjustment
  follows `mkr_adjust_project_id` and is omitted if only one project remains.

**Outputs and scope**

Results are under `output/{join.name}_join_zoom/{zoom.name}/`. The cell list,
selection summary and requested settings are saved as `selected_cells.csv.gz`,
`selection_summary.csv` and `effective_config.yaml`. Metadata and parent labels
are retained for comparison: the new Harmony PCs use `hmny_pca_*`, while parent
PCs use `parent_hmny_pca_*`. Other parent clusters and coordinates also have a
`parent_` prefix; `parent_label` records the selection label.

The report is `public/{join.name}_join_zoom_{zoom.name}.html`; the app defaults
to `public/shiny_{join.name}_join_zoom_{zoom.name}`. Cluster annotations must be
reviewed for the new clusters: provide a new `shiny.annotation_csv` when ready.
Any custom `shiny.app_tag` must differ from the parent and other zooms.

Switching marker resolution to another already calculated resolution reuses
selection, HVGs, PCA and integration. Changing selected cells or QC filters
requires recalculating the subset. Check the dry run when resuming an older run.

Optional cell-cycle correction is available in joins and joined zooms using
`integration.int_pca_method: bpcells`. Set `int_cell_cycle_correction` to
`within_sample` or `global`; it is omitted by default and inherited from the
parent join unless overridden. Every selected cell must have a finite
`tricycle_theta` score from its source project. Review phase calibration across
projects before correction. The subset fits its own correction model; marker
counts remain unchanged.

Joined zooms currently do not support recursive zooms, pool exclusions,
XGBoost training, fresh cell-type labelling or per-batch
subset H5AD/SCE exports.

### Configuration file

The join configuration reuses the `hvg`, `integration`, `marker_genes`, `label_celltypes`,
`train_xgboost`, `shiny`, and `resources` sections from a project configuration where supported by
the join schema. Use the current generated template and schema for exact fields and defaults. The
join must specify a `ref_txome`, or a `probe_set` for Flex data, that matches every source project;
only compatible projects can be joined.

Join integration may use arrays for `int_batch_var` and `int_theta` to configure
multiple Harmony covariates. When both values are arrays, they must have the same
length so that each batch variable has a corresponding theta value.

Additional parameters include:

* `name`: short name under `join`; output directories include `{name}_join` and `{name}_marker_genes`.
* `metadata_vars` (optional): list of column names from the source projects' `sample_metadata` CSV
  files to carry through into the joint sample metadata and presentation outputs. Each variable must
  exist in at least one project's metadata file, but does not need to be present in all projects.
  Missing values are filled with `NA`.
* `projects`: mapping containing at least two source projects. Each entry requires `config`, the
  path to a completed project config, and may specify `zoom_name` to join a completed zoom instead
  of the full project.
* `int_pca_method`: `bpcells` (default) reads normalized and scaled counts from
  disk for PCA, using a compact BPCells matrix and node-local storage. This avoids
  holding the full expression matrix in memory. `scanpy` calculates PCA in memory
  on CPU or GPU.

Join markers use the [same edgeR.bp method](#marker_genes) as projects and zooms.
Raw counts are summed by sample and cluster into a BPCells pseudobulk matrix.
Preparation filters genes and pseudobulks and calculates TMM normalization
factors. Marker testing and variable-gene ranking can then run independently.
Use `marker_genes.mkr_adjust_project_id: true` to adjust for project differences
in dispersion estimation and marker comparisons; this takes longer than the
default unadjusted analysis.

`resources.n_join_marker_genes_calc` controls marker-fitting CPUs; default `8`.
Each worker reads the disk-backed counts rather than receiving a full matrix
copy. More workers can reduce runtime but require more memory. Completed fits
are saved as `.fitted.rds` checkpoints, allowing a failed summary step to resume
without refitting when the inputs, settings and fitting code are unchanged.
Changed inputs or an unreadable checkpoint require a fresh fit.

The workflow provides the BPCells and `edger.bp` R packages. Creating the marker
environment for the first time requires access to `raw.githubusercontent.com`
for the [edger.bp Conda packages](https://github.com/marusakod/edger.bp/tree/conda-channel).

Joint HVGs are selected from the final HVG lists of the source projects. Genes
are ranked first by the number of projects in which they were selected as an
HVG, then by their median final rank among those projects. This preserves each
source project's sample- or group-aware HVG calculation while giving every
project equal weight.

When a metadata column is numeric in one source project and text in another, the
joint column is stored as text. This supports identifier columns such as
`patient_id` that use numeric IDs in one study and prefixed IDs in another.

The join HTML report includes distributions of post-filter cell-level QC metrics
(UMI count, detected genes, mitochondrial proportion, and spliced proportion) for
clusters at the `marker_genes:mkr_sel_res` resolution. QC values are collected from
the source projects and restricted to the cells included in the joint analysis.

Source sample metadata is combined by a separate lightweight join rule. If only
descriptive columns such as diagnosis, tissue, or therapy change, rerunning the join
refreshes the joint sample metadata and HTML report without rerunning the joint count
matrix, PCA, integration, clustering, marker genes, or cell-type labels. This does
not apply when changing identifiers or variables used analytically (for example a
batch variable or `hvg_metadata_split_var`); those changes require the affected
analysis steps to be rerun.

Example {{scjoin}} configuration file:

```yaml
join:
  name: my_join
  proj_dir: /path/to/join_output
  date_stamp: "2025-01-15"
  ref_txome: human_2024

  # For Flex, replace ref_txome with probe_set:
  # probe_set: human_v1
  # tenx_assay_type: flex
  your_name: Testy McUser
  affiliation: where you work
  metadata_vars: [Condition, Sex, Region]
projects:
  project_a:
    config: /path/to/config_a.yaml
    zoom_name: T_cells
  project_b:
    config: /path/to/config_b.yaml
    zoom_name: T_cells
  project_c:
    config: /path/to/config_c.yaml
hvg:
  hvg_n_hvgs: 2000
integration:
  int_pca_method: bpcells
  int_embedding: harmony
label_celltypes:
  - labeller: celltypist
    model: Immune_All_Low
train_xgboost:
  annots_f: /path/to/annotations.csv.gz
  xgboost_tag: my_classifier
shiny:
  app_tag: my_join
  app_title: My Joint Analysis
```

**Notes on `label_celltypes` in join:**

* When `label_celltypes` is configured, it runs as part of the default `scprocess join` target.
* Source projects' naive predictions can be reused for the same `labeller` and `model`. For XGBoost,
  the saved predictions must also match the installed model checksum; missing or different model
  checksums require fresh predictions. Labels are then aggregated using the join clusters.
* `save_cluster_names_file: true` generates a `cluster_names_for_shiny_*.csv` at the
  `marker_genes:mkr_sel_res` resolution. This file can be used as `annotation_csv` in the `shiny:`
  section.

The join report's **QC metrics by cluster** includes source scDblFinder scores,
shown on a logit scale with labels in the original score units. Scores are
matched by source project and cell ID using each project's merged real-cell
diagnostics. Missing diagnostic inputs or unmatched scores stop report
preparation with an error. Rendering the report does not rerun doublet detection.
