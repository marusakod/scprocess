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

Only one FASTQ source is allowed: use either `fastq_dir`, or `arv_uuids` with
`arv_instance`. For Flex, replace `ref_txome` with `probe_set` and set
`tenx_assay_type: flex`. Omit `train_xgboost` unless you intend to train a
classifier; run that analysis with `-r train_xgboost`. The example keeps the
original training labels: cluster majority voting is disabled unless explicitly
requested. For independent testing, reserve test groups as described under
[train_xgboost](#train_xgboost).

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

Quality control removes cells with measurements that are unsuitable for the
analysis. Common warning signs include very few detected RNA molecules or genes,
or a high mitochondrial proportion. These measurements also vary between real
cell types and between cells and nuclei, so review their distributions before
choosing cutoffs. A low-RNA population should not be rejected solely because
it differs from the most abundant cells.

This stage also identifies **doublets**: droplets containing material from more
than one cell. Individual doublet calls are made here; a later integration step
can remove whole clusters enriched for these calls. Review both stages when
assessing how many cells are retained.

###### Filtering cells and samples

UMIs count detected RNA molecules. Detected features usually correspond to
genes. Minimum count and feature thresholds remove cells with little usable
expression information; optional upper limits can flag unusually large
libraries, although a high count alone does not establish that a cell is a doublet.

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

Review the QC report's distributions and retained-cell counts by sample.
Check whether a proposed cutoff removes a plausible cell population or affects
one sample disproportionately. Count and feature limits include their boundary
values; a configured maximum must be at least as large as the minimum.

??? note "Which genes contribute to the measurements"

    * `exclude_mito`: exclude mitochondrial genes from the expression matrix
      used for subsequent analysis; default `true`. Mitochondrial content is
      still measured for QC.
    * The mitochondrial and splicing summaries are calculated from molecule
      counts. Small count adjustments keep proportions away from exactly 0
      and 1 for plotting and filtering on a logit scale.

###### Identifying individual doublets

`scDblFinder` compares observed cells with artificial doublets made by combining
expression profiles. It assigns each cell a score and a singlet/doublet call.
The calling threshold can differ between samples, so compare the scores,
thresholds and call rates together. Scores and call rates alone do not establish
whether a rare population is a genuine cell type or a doublet population.

The QC report shows observed and artificial-doublet score distributions, with
linear- and log-count histogram tabs and up to 20 samples per page. The
integration report shows where doublet-enriched clusters lie among the other
cells. Review both before changing the settings below.

??? note "Changing scDblFinder settings"

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

    The doublet-input filters act before calling doublets and are separate from
    the final `qc_min_counts` and other QC thresholds. Changing them changes the
    cells available to scDblFinder and can change its inferred doublet rate.
    Final count and feature maxima are applied after scDblFinder has run.
    Input filters, calling settings and any overrides are saved in the doublet
    summaries. Changed calls can affect HVG selection and integration too.

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

    Unlisted samples inherit the project setting. Explicit `null` uses the
    installed scDblFinder default, even if the project setting is numeric.
    For HTO/custom multiplexing, doublet calling is performed per pool: use
    `pool_id` instead. Sample-level `dbl_dbr_sd` overrides in these pooled
    workflows are rejected.


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

This stage groups cells with similar expression profiles and makes the UMAP
used to explore the dataset. It first summarizes variation across selected
genes using principal component analysis (**PCA**). It can then reduce
differences associated with sample preparation or processing using **Harmony**,
before finding clusters and calculating the UMAP.

Clustering and the UMAP serve different purposes: clusters group cells using
their expression relationships, while the UMAP displays those relationships
in two dimensions. A visually separate island is not by itself evidence of a
distinct cell type. Check marker expression and sample composition alongside
the plot.

###### Choosing which differences to correct

* `int_embedding`: use `harmony` (default) to correct for the chosen batch
  variable, or `pca` to cluster without batch correction.
* `int_batch_var`: the variable identifying batches for Harmony; default
  `sample_id`. For HTO/custom multiplexing, `pool_id` is also an option.
  Join analyses can correct for multiple variables, as described in the
  [join configuration](#configuration-file_2).
* `int_theta`: controls how strongly Harmony favours mixing batches;
  default `0.1`. Larger values encourage more mixing.

Choose batch variables according to the experimental design. For example,
sample preparation can introduce differences between libraries, but disease,
region and treatment may also differ between those same libraries. Inspect
known cell populations and their markers to check that correction has not
removed the differences you want to study. Improved sample mixing alone does
not establish that the correction is appropriate.

###### Choosing clusters and reviewing the result

* `int_n_dims`: number of principal components used to describe expression
  differences; default `50`.
* `int_res_ls`: clustering resolutions to calculate; default
  `[0.1, 0.2, 0.5, 1, 2]`. Larger resolutions generally produce more, smaller
  clusters. Calculating several lets you compare broad cell classes with
  finer subdivisions.

The resolution used for marker-gene analysis is selected separately with
`marker_genes.mkr_sel_res`. Review the UMAP coloured by sample, cluster, known
markers and available cell-type labels. Look for expected populations,
sample-specific groups, and clusters with poor QC or high doublet scores.
Investigate unexpected groups before deciding whether to keep or remove them.

###### Removing doublet-enriched clusters

Some mixed-cell profiles remain difficult to identify one at a time. A
preliminary clustering therefore checks for groups containing a high proportion
of cells already labelled as doublets. Cells in clusters above the
chosen cutoff are excluded before calculating the final integration.

* `int_dbl_res`: resolution of this preliminary clustering; default `4`.
  It is separate from the resolutions used to examine the final cell populations.
* `int_dbl_threshold_method`: `adaptive` (default) estimates the cutoff from
  the observed cluster doublet proportions. `manual` uses your chosen cutoff.
* `int_dbl_cl_prop`: cutoff for the manual method, expressed as a fraction.
  For example, `0.3` excludes clusters with **more than 30%** called doublets.
  It is required with `manual`, has no default, and must be omitted with
  `adaptive`.

For a manual cutoff:

```yaml
integration:
  int_dbl_threshold_method: manual
  int_dbl_cl_prop: 0.3
```

Review the integration report's doublet proportions and cluster QC plots.
Whole-cluster filtering can remove cells that were individually called
singlets, so inspect the biological identity and size of affected clusters.
It does not change the original scDblFinder calls. This cutoff is separate
from the per-cell calling threshold and `qc.dbl_dbr_sd`.

The applied cutoff and method are saved during integration and displayed in
the report. Changing them requires rerunning preliminary classification and
downstream integration; editing the configuration does not change saved results.
Reports require the threshold metadata from the integration run.

??? note "How the adaptive cluster cutoff is estimated"

    Clusters below 50% doublets and clusters above 50% define two groups.
    Clusters exactly at 50% are omitted from estimation. Each cluster contributes
    equally, regardless of its cell count. On the corrected logit scale, the
    median and scaled median absolute deviation (MAD; multiplier 1.4826) of each
    group define the centre and standard deviation of a Gaussian component.
    Mixing weights are the proportions of clusters in the two groups.

    The cutoff is the crossing of the two weighted densities between their
    centres, converted back to a doublet proportion. There is no 50% cap.
    Estimation stops with an error if either group contains fewer than two
    clusters, has zero MAD, or lacks a unique crossing between the centres.
    Review the cluster proportions and choose a manual cutoff if the adaptive
    estimate cannot be calculated or is unsuitable for the dataset.

??? note "UMAP layout settings"

    These settings change how cells are displayed and can usually be left at
    their defaults.

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

??? note "Computing options and saved objects"

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

    With BPCells, projects first calculate PCA using QC-passed cells together
    with observed cells called doublets (`pca_with_doublets`), then for retained cells alone
    (`pca_singlets_only`). Cell-cycle correction, if enabled, is fitted separately
    at each stage. Regression tables record its summaries and gene coefficients.


##### cell_cycle

Use this analysis to examine where cells lie along the cell cycle and whether
cell-cycle activity contributes to the structure of your dataset. It uses
`tricycle`, which compares expression of cell-cycle genes with a fixed reference.
Cycling cells can form a loop in the resulting plot. A cell's position around
this loop gives its estimated cell-cycle position, `tricycle_theta`.

Add the following block to estimate cell-cycle position and make the report:

```yaml
cell_cycle:
  cyc_density_percentile: 1.3
```

Run `scprocess run CONFIG -r cell_cycle`. The report is saved as
`public/{short_tag}_cell_cycle.html`. Leaving out the `cell_cycle` block disables
this analysis. Adding it does **not** remove cell-cycle effects from clustering.

###### Choosing and checking the centre of the loop

Cell-cycle position is measured around a centre, called the **origin**. Many
datasets contain a large group of non-cycling or G0/G1 cells and relatively few
cells in S or G2/M. To avoid that dense group dominating the estimated centre,
{{sc}} gives cells in less populated parts of the plot more weight.

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

Review the centre together with the expression of known cell-cycle markers.
The report shows how the centre changes across alternative settings, marker
expression around the loop, and results for each sample. A dataset with few
cycling cells may have no clear loop. Approximate phase labels and percentages
should therefore be interpreted alongside the marker plots, rather than as
definitive measurements of proliferation.

??? note "Additional plot and centre-estimation settings"

    These settings can usually be left at their defaults.

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

    The reference ellipse and the percentage outside it describe the spread of
    the cells. They do not determine the centre or identify cycling cells.
    Phase-percentage plots omit G0/G1 to make the less abundant phases easier
    to see.

###### Removing cell-cycle effects from clustering

Cell-cycle activity can separate otherwise similar cells into different
clusters. If this obscures the cell types or states you want to study, you can
remove the estimated cell-cycle contribution before PCA and clustering:

```yaml
cell_cycle:
  cyc_density_percentile: 1.3
integration:
  int_pca_method: bpcells
  int_cell_cycle_correction: within_sample
```

`int_cell_cycle_correction` has two options:

* `within_sample`: removes cell-cycle-associated differences between cells
  within each sample, preserving each sample's average expression on the scale
  used for PCA. Start here when you want to reduce cell-cycle structure before
  Harmony integration. Differences between samples with more or fewer cycling
  cells are retained at this step.
* `global`: also removes the predicted differences between sample averages
  attributable to their different cell-cycle composition. Consider this when
  those differences dominate the leading PCs. It does not remove all
  differences between samples.

Omit `int_cell_cycle_correction` to leave expression uncorrected. Correction
requires cell-cycle estimates and `int_pca_method: bpcells`. It changes the
expression used for PCA, clustering and UMAP, but leaves the count matrices
used for differential expression and the estimated cell-cycle positions
unchanged. Different clusters can nevertheless change which cells contribute
to each pseudobulk.

Compare corrected and uncorrected results using known markers, cell types,
sample identities and experimental conditions. Better sample mixing alone
does not establish a better result. If proliferation is part of the biological
question, correction could remove information you want to retain.

??? note "Correction settings and assumptions"

    * `int_cell_cycle_harmonics`: how flexibly the correction follows changes
      around the cycle; `1` or `2`, default `2`.
    * `int_cell_cycle_ridge_lambda`: limits the strength of the fitted
      correction; default `0.1`. Larger values reduce the correction.

    Both belong in `integration` and are accepted only when
    `int_cell_cycle_correction` is set. They can usually be left at their
    defaults.

    The model estimates a shared relationship between expression and cell-cycle
    position using variation within samples. It assumes that this relationship
    is sufficiently similar across samples. Other biological changes that
    correlate with cell-cycle position can also be affected. Both correction
    options use this same fit; they differ in whether predicted shifts between
    sample averages are removed. Harmony cannot recover information already
    lost during PCA, so the options need not give the same final clusters.

    In projects, correction is fitted separately before preliminary and final
    PCA. Project zooms fit it for their selected cells. Cells without a sample
    assignment in HTO/custom data use their physical run for centring.

??? note "Calculation details and saved outputs"

    The species is inferred from `project.ref_txome` or `project.probe_set`.
    Projections are calculated per physical run: `sample_id` for
    non-multiplexed, Flex and OCM data; `pool_id` for HTO/custom data. All runs
    use the same reference, and coordinates are centred once across the whole
    project. The centre is estimated from QC-passed candidate singlets; known
    doublets receive coordinates and positions but do not influence its estimate.

    The density threshold is the requested percentile of the eligible cells'
    density values, with linear interpolation between adjacent values. A cell's
    weight is `min(1, density_threshold / cell_density)`. The centre is the
    weighted mean of its coordinates. Percentile `0` uses the minimum density,
    giving weights proportional to inverse density; percentile `100` uses the
    maximum density, giving every cell weight 1. The saved summaries record the
    threshold and the actual percentage of cells at full weight. Diagnostic
    display samples use these weights as sampling probabilities, without
    changing the fitted centre.
    The reference ellipse uses minimum covariance determinant (MCD) and a 99%
    boundary; the percentage outside is a diagnostic, not a cycling-cell estimate.

    Saved results include:

    * `tricycle_pc1` and `tricycle_pc2`: coordinates in the fixed cell-cycle
      reference, after project-wide centring. These are different from the PCs
      used for integration.
    * `tricycle_theta`: position around the chosen centre in radians, from
      `0` up to but excluding `2*pi` (equivalent to 0–360 degrees).

    Coordinates and positions depend on project composition. Changing the
    weighting or manual centre reuses the per-run projections, but requires
    recalculating positions and dependent results.

    Results are in `output/{short_tag}_cell_cycle/`; the report source is
    `analysis/{short_tag}_cell_cycle.Rmd`. The report includes expression of
    12 cell-cycle marker genes, read from the full count matrices independently
    of HVG selection. Spliced, unspliced and ambiguous features for the same
    symbol are summed before normalization. Missing genes are recorded in the
    run summary and retained as zero-valued columns in
    `tricycle_marker_expression_*.csv.gz`.

    Correction summaries and gene coefficients are saved when correction is
    enabled. Predicted changes in sample averages should be approximately zero
    for `within_sample`. If cell-cycle composition needs to be accounted for
    in a separate differential-expression analysis, use a model appropriate to
    that comparison; circular positions require both sine and cosine terms,
    rather than an ordinary average angle or mean sine alone.

??? note "Cell-cycle resource settings"

    Per-run projection uses `calculate_tricycle_projection_one_run`. Override
    memory and runtime with `gb_calculate_tricycle_projection_one_run` and
    `mins_calculate_tricycle_projection_one_run`.

    Pooling and the reference ellipse each request 8 GB and 30 minutes by
    default. Their resource keys are `gb_aggregate_tricycle`,
    `mins_aggregate_tricycle`, `gb_estimate_tricycle_reference_ellipse` and
    `mins_estimate_tricycle_reference_ellipse`.


##### marker_genes

Marker analysis asks **which genes distinguish each cluster from the other
clusters**. Use these genes to help identify cell populations and investigate
their biology. For example, several consistent astrocyte markers provide
stronger support for an astrocyte annotation than a single highly expressed gene.

For each cluster, {{sc}} sums counts separately within each sample. These
sample-level expression profiles are called **pseudobulks**. The marker test
compares each cluster with all remaining clusters, accounting for variation
between samples. Projects, joins and zooms use the same edgeR.bp calculation.
This comparison identifies cluster-associated genes; a treatment or disease
comparison requires a separate differential-expression analysis with a model
appropriate to the experimental design.

###### Selecting the clusters and minimum cell counts

* `mkr_sel_res`: clustering resolution used for marker analysis; default `0.2`.
  Choose one of the resolutions already calculated in `integration.int_res_ls`.
* `mkr_min_cl_size`: minimum total cells in a cluster to calculate its markers;
  default `100`.
* `mkr_min_cells`: minimum cells in each cluster–sample pseudobulk for it to
  contribute to marker testing; default `10`.

These minimums apply to marker analysis, rather than removing those cells from
the integration. Small populations may have few eligible samples after
filtering; inspect the sample contributions when interpreting their markers.

###### Interpreting marker results

A positive **log fold change** (`logFC`) means higher expression in the selected
cluster than in the comparison group; a negative value means lower expression.
The **false discovery rate** (`FDR`) is the P value adjusted for testing many
genes. Inspect the size of the expression difference, its statistical support
and the expression levels across samples together. A very small FDR does not
necessarily make a gene a useful or specific biological marker.

The report shows marker expression, sample-level pseudobulks and variable
genes. Check whether markers agree with one another and whether the apparent
difference is supported across samples rather than dominated by a single sample.
The comparison is against the other clusters present in this analysis, so a
marker list can change when the selected cell populations or clustering change.

* `mkr_min_cpm_mkr`: minimum expression for selecting genes shown in marker
  plots; default `50` counts per million (CPM). This is an expression filter
  for the report, not a significance cutoff for the marker test.
* `mkr_not_ok_re`: pattern matching gene types to omit from marker-expression
  plots; default `"(lincRNA|lncRNA|pseudogene|antisense)"`.

###### Gene-set enrichment

Gene-set enrichment analysis (**GSEA**) asks whether groups of related genes
are concentrated towards the higher- or lower-expression end of a cluster's
ranked marker results. It can help summarize biological processes, but an
enriched gene set does not by itself establish that a pathway is activated.

* `mkr_do_gsea`: run GSEA; default `true`.
* `mkr_gsea_cut`: FDR cutoff used for GSEA results; default `0.1`.

??? note "Additional GSEA settings"

    * `mkr_min_cpm_go`: minimum expression for genes included in GSEA;
      default `1` CPM.
    * `mkr_max_zero_p`: maximum proportion of a cluster's pseudobulks with
      zero counts when filtering genes for enrichment; default `0.5`. The zero
      filter is not applied to genes with negative `logFC`.
    * `mkr_gsea_var`: how genes are ranked. The default, `z_score`, uses a
      z-score based on signed `log10(FDR)`, combining the direction and
      statistical support of the difference. `logFC` ranks by fold change.

??? note "Adding custom marker-gene panels"

    `mkr_custom_genesets` supplies named gene panels for the report. Each entry
    contains:

    * `name`: a unique name for the panel.
    * `file`: path to a CSV containing `label` (the marker category), and
      `symbol` and/or `ensembl_id`. If the path is omitted, {{sc}} looks for
      `$SCPROCESS_DATA_DIR/marker_genes/{name}.csv`.

??? note "Normalization and variation between samples"

    Counts are normalized using TMM to account for differences in library size
    and RNA composition. Each cluster is compared with the remaining clusters
    using an edgeR Treat test.

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

    Variable genes in marker reports are ranked by variance of TMM-normalized
    logCPM across eligible cluster–sample pseudobulks. Profiles containing only
    one cell are excluded from this ranking; marker testing uses the separate
    `mkr_min_cells` criterion.

??? note "Files, resource settings and changes from earlier versions"

    Pseudobulks are saved as `pb_*.bpcells` directories and read in blocks,
    avoiding the need to load the full count matrix into memory. Preparation,
    marker fitting and plot data are saved separately, so a failed stage can
    resume without rebuilding the counts.

    `resources.n_run_marker_genes` controls fitting CPUs for projects and
    project zooms. Fitting memory and time use `gb_run_marker_genes` and
    `mins_run_marker_genes`, or the corresponding `zoom_` rule names. Other
    marker stages have separate [resource settings](#resources).

    Since 0.5.0, projects and zooms use the same normalization, dispersion and
    variability calculations as joins. Their pseudobulks are BPCells directories
    rather than RDS objects; marker CSV names are unchanged. These methodological
    changes can alter results relative to earlier project-marker calculations.
    The former combined in-memory resource model is no longer used.


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
* `label_map`: path to a CSV file mapping fine-grained labels to coarse labels. Must contain columns
  `label` and `coarse_label`. Should be absolute or relative to `proj_dir`. If not specified, {{sc}}
  will use the default label map included with the classifier (if available).
* `save_cluster_names_file`: if `true`, generate a `cluster_names_for_shiny_*.csv` file mapping
  clusters (at the `marker_genes:mkr_sel_res` resolution) to predicted cell type names. This file
  can be used as the `annotation_csv` in the shiny config. Requires a `marker_genes` block to be
  configured. Default is `false`.


##### train_xgboost

Use this analysis to train a classifier from cells you have already labelled,
so that the same labels can be assigned to other datasets. For example, a
classifier trained on a labelled brain atlas can predict cell types in a new
brain study. Its usefulness depends on the quality of those labels and how
well the training cells represent the cells you later want to classify.

Training uses expression and integration results from a completed project,
join or project zoom. Supply a table linking cell IDs to labels, then add a
`train_xgboost` block. For example, if your metadata includes `donor_id`:

```yaml
train_xgboost:
  annots_f: data/configs/cell-labels.csv
  validation_group_var: donor_id
  validation_fraction: 0.2
  train_cells_per_label: 2000
  validation_cells_per_label: 500
```

Run `scprocess run CONFIG -r train_xgboost`, or
`scprocess join CONFIG -r train_xgboost` for a join. Review the resulting HTML
report before installing the classifier for use in other projects.

###### Input labels

* `annots_f` (required): CSV or CSV.gz with a unique `cell_id` for each labelled
  cell and a `label` column, such as `Astrocyte` or `Microglia`. The file can
  also supply donor, region or other metadata used to choose cells. Values
  shared with the integration metadata must agree. Relative paths are resolved
  from the project folder.
* `excluded_labels`: labels to leave out of training and the main evaluation;
  default `[]`. For example, you might exclude an uncertain annotation rather
  than teach the classifier that it is a distinct cell type.
* `min_cells_per_label`: minimum number of eligible cells needed to retain a
  label before sampling; default `20`. Labels with fewer cells are excluded.
  This minimum does not by itself establish that a rare label has enough
  examples for reliable classification.

Original labels are kept unless you explicitly enable cluster majority voting,
described below.

###### Choosing cells for training and checking predictions

Cells have different roles in the analysis:

* **Training cells** teach the classifier the expression patterns associated
  with each label.
* **Validation cells** come from separate groups, such as donors or tissue
  tubes. Their predictions are checked during training to decide when to stop.
  Because these cells influence that decision, agreement with their labels can
  overestimate performance on a new dataset.
* **Test cells**, when requested, come from groups held out from both steps.
  They provide a more independent check of predictions.
* **Remaining cells** were not selected because of the limits on training and
  validation cell numbers. They can share donors or samples with training
  cells, so they are not an independent test.

Choose groups that match the intended use. To assess predictions in new donors,
keep whole donors together rather than separating cells from the same donor
between training and validation.

* `validation_group_var`: metadata column defining the groups to keep together,
  for example `donor_id` or a tissue-tube identifier. By default it uses the
  processing batch inherited from the project or join.
* `validation_fraction`: fraction of groups assigned to validation; default
  `0.2`. With ten donors, this initially selects two donors, not 20% of cells.
  The final cell fraction depends on how many cells each donor contributes.
* `test_group_var` and `test_groups`: column and values to hold out for testing.
  For example, `test_group_var: donor_id` and `test_groups: [donor_3]` reserve
  that donor. Supply both settings or neither. Test groups are selected first.
* `train_cells_per_label` and `validation_cells_per_label`: maximum cells to
  select per label after assigning groups; defaults `2000` and `500`. These
  limits keep abundant labels from dominating simply through their cell counts.

The report compares predicted and supplied labels, shows their distribution
over the UMAP, and summarizes how often each label is recovered or confused
with another. Examine rare labels and biologically plausible confusions,
alongside the overall accuracy. Agreement on training cells is expected to
be higher than agreement on genuinely new samples.

??? note "Balancing regions, donors or other metadata"

    * `validation_variable_combinations`: choose validation groups separately
      within combinations of metadata values. For example, `[region, condition]`
      can help represent each region/condition combination on both sides.
      Every group must have one consistent combination of these values, and
      enough groups must be available to split. Omit for a single overall split.
    * `sampling_variable_combinations`: balance cell selection within each label
      across combinations such as `[donor_id, region]`. This happens after group
      assignment and separately for training and validation. Omit to give cells
      within each label and assigned split an equal chance of selection.

    `validation_fraction` is rounded down to a whole number of groups, with at
    least one on each side where two or more groups are available. A combination
    with only one group contributes none to validation. Whole groups can be
    moved back to training to retain enough cells for rare labels, so inspect
    the saved split counts rather than assuming the requested fraction was met.
    Training stops with an error if no suitable validation groups remain.

###### Label names, colours and coarse groups

Two optional CSV files control how labels appear in reports:

* `label_specs_f`: one row per label, in the desired display order. Required
  column: `label`. Optional columns: `colour` and `coarse_label`.
* `coarse_label_specs_f`: one row per coarse group, in its display order.
  Required column: `coarse_label`; optional column: `coarse_colour`. This file
  requires a coarse mapping in `label_specs_f`.

For example, `Astrocyte` and `Bergmann glia` can both map to `Astroglia`.
The classifier still predicts the original labels; the report can also compare
the corresponding coarse groups. Each coarse group is a direct lookup from
the predicted label. It does not combine probabilities or vote across cells.

Include every model label in `label_specs_f`. If colours are supplied, give
every label a hex colour code such as `#74C476`. Omit coarse mappings to show
fine labels only, or omit both files to use the default display. Changing
order, colours or coarse mappings updates the report without retraining the
model or recalculating predictions.

??? note "Replacing training labels by cluster majority vote"

    This option can reduce inconsistent labels within a cluster, but can also
    overwrite genuine rare cell types that share a cluster with more abundant
    cells. Leave it off when preserving those distinctions is important.

    * `cluster_majority_vote_resolution`: an already calculated clustering
      resolution at which to take each cluster's most frequent training label.
      Omit to keep the original labels.
    * `cluster_majority_vote_min_purity`: minimum fraction of labelled training
      cells that must agree before their label can be assigned to the cluster;
      default `0.9`. For example, at least 90% must share the most frequent label.

    At least ten training annotations are also needed per cluster. This is a
    fixed minimum, not a configurable setting. Votes are determined from
    training cells, keeping validation and test labels out of that decision.

??? note "Gene selection and training settings"

    The classifier is trained twice. The first fit considers all genes remaining
    after filtering and identifies those most useful for distinguishing labels.
    The second fits the final model using the selected genes. A gene's **gain**
    measures its contribution to improving predictions during training; it is
    not a statistical test of differential expression.

    * `min_cells_expressed`: a gene must be detected in at least this many
      training cells; default `10`.
    * `gene_exclude_re`: pattern matching gene types to exclude; default
      `"(lincRNA|lncRNA|pseudogene|antisense)"`. Excluding these features can help
      when applying a model to assays that do not measure them, but does not
      guarantee coverage of all selected genes in a new assay.
    * `gain_threshold`: select genes accounting for this fraction of total
      gain in the first fit; default `0.9` (90%).
    * `min_genes`, `max_genes`: lower and upper limits on the selected gene
      count; defaults `100` and `3000`.
    * `seed`: random seed for reproducible cell selection and fitting;
      default `42`.
    * `use_gpu`: use a GPU for training; default `false`.
    * `xgboost_tag`: name used in model filenames. Supplied names receive the
      `xgboost_` prefix. If omitted, the prefix is `xgboost_{full_tag}` for a
      project or `xgboost_{join.name}` for a join.

    XGBoost builds predictions by successively adding decision trees. The
    settings below control how it does this; most users can start with the
    defaults. Validation predictions determine when to stop adding trees.

    | Setting | Meaning | Default |
    |---------|---------|---------|
    | `pass1_subsample` | Fraction of training cells sampled for each round in the first fit | `0.632` |
    | `pass1_colsample_bytree` | Fraction of genes considered for each tree in the first fit | `0.1` |
    | `pass1_learning_rate` | How much each added tree changes predictions | `0.1` |
    | `pass1_nrounds` | Maximum rounds in the first fit | `300` |
    | `pass1_early_stopping` | Stop after this many rounds without improved validation performance | `10` |
    | `pass2_colsample_bytree` | Fraction of selected genes considered for each tree in the final fit | `0.5` |
    | `pass2_learning_rate` | How much each added tree changes predictions | `0.05` |
    | `pass2_nrounds` | Maximum rounds in the final fit | `500` |
    | `pass2_early_stopping` | Stop after this many rounds without improved validation performance | `10` |

Once reviewed, the classifier can be
[installed for use in other projects](#setup-xgboost).


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

Building a Shiny app has two stages. `build_shiny_data` (or
`build_zoom_shiny_data`) performs the expensive H5AD reads, UMAP subsampling,
BPCells writes, pseudobulking, and marker/GSEA preparation. The lightweight
`configure_shiny_app` (or `configure_zoom_shiny_app`) copies the application
template and writes `annotation.csv`, `home.md`, and `shinyconfig.yaml`. Running
`scprocess ... -r shiny` targets the configuration stage and brings in the data
stage only when its analytical inputs are stale. Changes to cluster names,
order or colours, `app_title`, `email`, `keyword`, `default_gene`,
`metadata_labels`, `metadata_combns`, `home_md`, `cluster_palette`, or
`metadata_palettes` therefore do not rebuild the expression matrices.

Changing `n_keep`, `metadata_vars`, clustering/marker resolution, sample
metadata, integration outputs, H5ADs, markers, HVGs, or GSEA inputs rebuilds the
Shiny data. Check the dry run before rebuilding a large app.

The cluster UMAPs include a **Repel cluster labels** option. Repelled label
positions are computed from the full UMAP during the app build: labels are moved
into low-density space, given clearance from cells, and allocated globally so
their numbered circles do not overlap one another. Rebuild the Shiny target to
refresh these positions after changing the integration or clustering.

The sample metadata CSV is an explicit input to the Shiny data stage. Editing
descriptive values therefore marks the generated Shiny data for rebuilding. For a join, first
run the normal join target so its joint sample metadata snapshot is refreshed, then
run `scprocess join CONFIG -r shiny` to rebuild the app.


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
