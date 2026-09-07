# Reference

## {{scsetup}} { #scprocess-setup }

**Description**: Download all data required for {{sc}} and index reference transcriptomes and probe sets for `simpleaf`.

**Parameters**:

* `-c`/`--rangerurl`:  download link for Cell Ranger (v9.0.0 or higher) available on the [10x Genomics CellRanger download & installation page](https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions); only required when running the command for the first time.

### configuration file

The command requires a configuration file named `scprocess_setup.yaml` located in {{sc}} data directory (for instructions on how to set up the {{sc}} data directory see the [Getting started](setup.md#scprocess-data-directory-setup) section). In this file, the user can specify parameters that are used across all {{sc}} projects, such as HPC configuration and reference genomes that will be made available for {{sc}}. For example:

```yaml
user:
  profile_template: slurm_default 
  profile_name:     slurm_my_cluster
  your_name:      Testy McUser
  affiliation:    Unemployed
  int_use_gpu:    false
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

* `profile_template`: the name of the bundled HPC profile template to copy from the {{sc}} _profiles_ directory into `$SCPROCESS_DATA_DIR/profiles`. Bundled templates include `slurm_default`, `lsf_default`, `pbs_default`, `sge_default`, and `htcondor_default`. Exactly one of `profile_template` and `local_cores` should be specified.
* `profile_name` (optional): name of the active profile directory under `$SCPROCESS_DATA_DIR/profiles`. Defaults to the value of `profile_template`.
* `profile` (deprecated): accepted as an alias for `profile_template` for existing setup files.
* `local_cores`: number of CPU cores available for local execution (see [Snakemake documentation](https://snakemake.readthedocs.io/en/v9.8.0/executing/cli.html) for more details). Exactly one of `profile_template` and `local_cores` should be specified.
* `your_name` (optional): author's name. If specified it will be used in the configuration file for new projects created with the `scprocess newproj -c` command.
* `affiliation` (optional): author's affiliation. If specified it will be used in the configuration file for new projects created with the `scprocess newproj -c` command.
* `int_use_gpu` (optional): whether to use GPU acceleration (`RAPIDS-singlecell`) for integration and clustering steps. If `false` the value will be used in the configuration file for new projects created with the `scprocess newproj -c` command.

##### arvados

* `arv_instance` (optional): the name of the default Arvados instance for the user. If specified it will be used in the configuration file for new projects created with the `scprocess newproj -c` command.

##### ref_txomes

Prebuilt human and mouse reference transcriptomes from 10x Genomics can be downloaded with {{scsetup}} by adding `tenx` to the `scprocess_setup.yaml` file. Valid values for names are `human_2024`, `mouse_2024`, `human_2020`, `mouse_2020`.

Names and specifications for custom references should be listed in the `custom` section of the `scprocess_setup.yaml` file. For each `custom` genome users have to provide the following parameters:

* `name`: name to be used for the reference
* `fasta`: path to FASTA file
* `gtf`: path to GTF file
* `mito_str`: regular expression used to identify genes in the mitochondial genome (example for mouse: `"^mt-"`)

Optional parameters for both `tenx` and `custom` references are:

* `decoys`: whether or not poison k-mer information should be inserted into the index. This parameter is optional. If not specified, it defaults to `true` for all genomes.

Optional paramater for `tenx` references is:

* `rrnas`: whether or not ribosomal RNAs should be included in the reference. If not specified it defaults to `true` for all `tenx` genomes.

!!! note "Impact of custom parameters for `tenx` genomes on `scsetup` runtime"

    When configuring `tenx` genomes with their default values, {{scsetup}} will download prebuilt indices for `simpleaf`. However, if the default parameters are modified (e.g., setting `rrnas` or `decoys` to `false`), {{scsetup}} will build the indices from scratch during execution, which will increase the runtime.


!!! info "More about decoys"
    {{sc}} utilizes `simpleaf`, a lightweight mapping approach that, by default, maps sequenced fragments exclusively to the transcriptome. However, this can lead to incorrect mapping of reads that arise from unannotated genomic loci to the transcriptome. To mitigate this issue, the `decoys` parameter for `ref_txomes` is set to `true`. This option allows `simpleaf` to identify genomic regions with sequences similar to those in transcribed regions (decoys), thereby reducing the likelihood of false mappings. We strongly recommend keeping the decoy setting enabled. For further details, refer to Srivastava et al., 2019[@Srivastava2020-jb].

##### probe_sets

Probe set indices for 10x Flex data are prepared with {{scsetup}} by adding a `probe_sets` section to the `scprocess_setup.yaml` file. Each entry under `tenx` specifies a named probe set to download and index. Valid values for `name` are `human_v1`, `human_v2`, `mouse_v1`, and `mouse_v2`.


##### xgboost { #setup-xgboost }

Custom `XGBoost` classifiers trained with {{sc}} (using the `train_xgboost` functionality) can be registered with {{scsetup}} so they are available for cell type annotation in any project. To register a user-trained classifier, add an `xgboost` section to the `scprocess_setup.yaml` file. Each entry specifies a classifier to install. There are two ways to reference a trained classifier:

 * use the `config` parameter to point to the project configuration file that was used to train the classifier. {{sc}} will automatically locate the training outputs.
 * use the `scr_dir` parameter to point directly to the directory containing the trained classifier. Required the `name` parameter (unique identifier for the classifier) to be specified.

## {{scnew}}

**Description**: Create a new `workflowr` project directory for {{sc}} outputs.

**Parameters**:

* `name` (positional): name of the new `workflowr` project directory.
* `-w`/`--where` (optional): path to the directory where the new project will be created; defaults to the current working directory.
* `-s`/`--sub` (optional): if provided, creates `data/fastqs` and `data/metadata` subdirectories within the project.
* `-c`/`--config` (optional): generates a template configuration YAML file. If provided, it must be followed by either `sc` (single-cell) or `sn` (single-nucleus) to define standard QC thresholds. The generated config always includes `tenx_assay_type` in the `project` section. Additional modifiers can be appended:
    + `flex`: sets `tenx_assay_type: flex` and includes `probe_set` instead of `ref_txome`. e.g. `scprocess newproj project_name -c sc flex`
    + `multiplex`: adds a `multiplexing` section. When combined with `flex`, sets `demux_type: flex`; otherwise leaves `demux_type` blank for the user to fill in (e.g. `hto` or `custom`). e.g. `scprocess newproj project_name -c sc multiplex` or `scprocess newproj project_name -c sc flex multiplex`

## {{scnewjoin}} { #scprocess-newjoin }

**Description**: Create a new `workflowr` project directory for a {{scjoin}} outputs.

**Usage**:

```
scprocess newjoin <name> [-w <where>] [-p <config>...]
```

**Parameters**:

* `name` (positional): name of the new `workflowr` project directory.
* `-w`/`--where` (optional): path to the parent directory where the new project will be created; defaults to the current working directory.
* `-p`/`--projects` (optional): one or more paths to existing {{sc}} project configuration files.


## {{scknee}} { #scprocess-plotknee }

**Description**: Create an interactive barcode-rank plot. Can only be used once the mapping step is completed.

**Parameters**:

* `sample`: sample_id corresponding to the barcode-rank curve.
* `-k`/`--kneefile`: path to the knee plot data file generated by {{sc}}, e.g. `output/[short_tag]_mapping/af_[sample_id]/knee_plot_data_[sample_id]_[date_stamp].csv.gz`. Exactly one of `--kneefile` and `--configfile` should be specified.
* `-c`/`--configfile`: path to configuration file used for running {{sc}}. Exactly one of `--kneefile` and `--configfile` should be specified.


## {{scrun}} { #scprocess-run }

**Description**: Run {{sc}}.

**Parameters**:

* `configfile` (positional): path to a configuration YAML file.
* `-n`/`--dry-run` (optional): perform a trial run which lists all steps that {{sc}} would do and does not create any new files. Helpful for checking input files and parameters.
* `--create-envs` (optional): only create the conda environments needed for the workflow, without running any rules.
* `-E`/`--extraagrs` (optional): list of additional arguments to pass to `Snakemake`. Refer to [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for a detailed explanation of available command-line options.
* `-r`/`--rule` (optional): Specifies which rule {{sc}} should run. The options are:
    + `all`: default; includes all [Core pipeline steps](introduction.md#core-pipeline-steps) plus `label_celltypes` (if configured).
    + `mapping`: read alignment and quantification.
    + `ambient`: ambient RNA removal (optional) and cell calling.
    + `demux`: sample demultiplexing.
    + `qc`: qc filtering.
    + `hvg`: calculation of highly variable genes.
    + `integration`: dimentionality reduction with PCA, optional batch correction with `Harmony`, UMAP and clustering.
    + `marker_genes`: marker gene identification and optional gene set enrichment analysis.
    + `label_celltypes`: cell type annotation using a pre-trained classifier.
    + `train_xgboost`: train a custom `XGBoost` classifier from annotated data.
    + `zoom`: subclustering.
    + `shiny`: build an interactive Shiny app from pipeline outputs.
* `--zoom` (optional): only used with `-r shiny`. Build the Shiny app for a zoom (subclustering) output. Pass a zoom name (e.g. `--zoom immune_cells`) or `all` to build apps for every zoom defined in the config.


### configuration file

This is an example `configfile` for {{sc}} with all parameters and their default values/placeholders. Required parameters are highlighted:

=== "default values"

    ```yaml hl_lines="1 2 3 4 5 6 7 8 9 10 11 12 13 14"
    project:
      proj_dir:
      fastq_dir: # should not be defined if arv_uuids is defined
      arv_uuids: # should not be defined if fastq_dir is defined
      arv_instance: # should only be defined if arv_uuids is defined
      full_tag:
      short_tag:
      your_name:
      affiliation:
      date_stamp:
      sample_metadata:
      tenx_assay_type:
      ref_txome:
      probe_set:
      metadata_vars:
      show_arv_uuids: true
      custom_sample_params:
      tenx_chemistry:
      exclude:
        sample_id:
        pool_id:
    multiplexing:
      demux_type: none
      fastq_dir:
      arv_uuids:
      feature_ref:
      demux_output:
    ambient:
      ambient_method: decontx_background
      cell_calling: barcodeRanks
      cb_version: v0.3.2
      cb_max_prop_kept: 0.9
      cb_learning_rate:
      cb_posterior_batch_size: 128
      cb_empty_training_fraction:
      cb_expected_cells:
      cb_total_droplets_included:
      cb_low_count_threshold:
    qc:
      qc_min_counts: 500
      qc_max_counts:
      qc_min_feats: 300
      qc_max_feats:
      qc_min_mito: 0
      qc_max_mito: 0.1
      qc_min_splice: 0
      qc_max_splice: 1
      qc_min_cells: 100
      dbl_min_feats: 100
      exclude_mito: true
    pb_empties:
      ambient_genes_logfc_thr: 0
      ambient_genes_fdr_thr: 0.01
    hvg:
      hvg_method: sample
      hvg_n_hvgs: 2000
      hvg_exclude_ambient_genes: True
      hvg_exclude_from_file:
      hvg_chunk_size: 2000
      hvg_metadata_split_var:
    integration:
      int_pca_method: bpcells
      int_use_gpu: true
      int_embedding: harmony
      int_theta: 0.1
      int_batch_var: sample_id
      int_n_dims: 50
      int_dbl_res: 4
      int_dbl_cl_prop: 0.5
      int_sce_outs: false
      int_res_ls: [0.1, 0.2, 0.5, 1, 2]
      int_use_paga: true
      int_paga_cl_res: 2
    marker_genes:
      mkr_sel_res: 0.2
      mkr_min_cl_size: 100
      mkr_min_cells: 10
      mkr_not_ok_re: "(lincRNA|lncRNA|pseudogene|antisense)"
      mkr_min_cpm_mkr: 50
      mkr_min_cpm_go: 1
      mkr_max_zero_p: 0.5
      mkr_do_gsea: true
      mkr_gsea_cut: 0.1
      mkr_gsea_var: z_score
      mkr_custom_genesets:
      - name:
        file:
    label_celltypes:
      - labeller:
        model:
        hi_res_cl: "RNA_snn_res.2"
        min_cl_prop: 0.5
        min_cl_size: 100
    train_xgboost:
      annots_f:
      ref_tag:                            # optional; defaults to xgboost_{full_tag}
      label_map_f:
      refine_labels: true
      purity_threshold: 0.65
      n_cells_per_type: 1000
      min_cells_per_type: 20
      min_cells_expressed: 10
      gene_exclude_re: "(lincRNA|lncRNA|pseudogene|antisense)"
      seed: 42
      use_gpu: false
      pass1_subsample: 0.632
      pass1_colsample_bytree: 0.1
      pass1_learning_rate: 0.1
      pass1_nrounds: 300
      pass1_early_stopping: 10
      pass2_colsample_bytree: 0.5
      pass2_learning_rate: 0.05
      pass2_nrounds: 500
      pass2_early_stopping: 10
      gain_threshold: 0.9
      min_genes: 100
      max_genes: 3000
    zoom:
    shiny:
    resources:
      retries: 3
      n_run_mapping: 8
    ```

=== "placeholders"

    ```yaml hl_lines="1 2 3 4 5 6 7 8 9 10 11 12 13 14"
    project:
      proj_dir: /path/to/proj/directory 
      fastq_dir: /path/to/directory/with/fastq/files
      arv_uuids: ["arkau-qr8st-1a2b3c4d5e6f7g8", "arkau-9v0wx-h9i8j7k6l5m4n3o", "arkau-z2y3x-p0q1r2s3t4u5v6w"]
      arv_instance: instance_name
      full_tag: test_project
      short_tag: test
      your_name: Testy McUser
      affiliation: where you work
      date_stamp: "2050-01-01"
      sample_metadata: /path/to/metadata.csv
      tenx_assay_type: poly_a
      ref_txome: human_2024
      probe_set: human_v1
      tenx_chemistry: 3v3
      metadata_vars: [var1, var2, var3]
      show_arv_uuids: true
      custom_sample_params: /path/to/file/with/custom_parameters.yaml
      exclude:
        sample_id:
          - sample1
          - sample2
        pool_id:
          - pool1
          - pool2
    multiplexing:
      demux_type: hto
      fastq_dir: /path/to/directory/with/hto_fastq/files
      arv_uuids: ["arkau-qr8st-1a2b3c4d5e6f7g8", "arkau-9v0wx-h9i8j7k6l5m4n3o", "arkau-z2y3x-p0q1r2s3t4u5v6w"]
      feature_ref: /path/to/feature_ref.csv
      demux_output: /path/to/demux_output.csv
    ambient:
      ambient_method: cellbender
      cb_version: v0.3.2
      cb_empty_training_fraction:
      cb_expected_cells: 10000
      cb_total_droplets_included: 20000
      cb_low_count_threshold: 5
      cb_learning_rate: 0.001
      cb_posterior_batch_size: 128
    qc:
      qc_min_counts: 500
      qc_max_counts:
      qc_min_feats: 300
      qc_max_feats:
      qc_min_mito: 0
      qc_max_mito: 0.1
      qc_min_splice: 0
      qc_max_splice: 1
      qc_min_cells: 100
      dbl_min_feats: 100
      exclude_mito: true
    pb_empties:
      ambient_genes_logfc_thr: 0
      ambient_genes_fdr_thr: 0.01
    hvg:
      hvg_method: sample
      hvg_n_hvgs: 2000
      hvg_exclude_ambient_genes: True
      hvg_exclude_from_file: /path/to/file/with/genes/to/exclude
      hvg_chunk_size: 2000
      hvg_metadata_split_var: var1
    integration:
      int_pca_method: bpcells
      int_use_gpu: true
      int_embedding: harmony
      int_theta: 0.1
      int_batch_var: sample_id
      int_n_dims: 50
      int_dbl_res: 4
      int_dbl_cl_prop: 0.5
      int_sce_outs: false
      int_res_ls: [0.1, 0.2, 0.5, 1, 2]
      int_use_paga: true
      int_paga_cl_res: 2
    marker_genes:
      mkr_sel_res: 0.2
      mkr_min_cl_size: 100
      mkr_min_cells: 10
      mkr_not_ok_re: "(lincRNA|lncRNA|pseudogene|antisense)"
      mkr_min_cpm_mkr: 50
      mkr_min_cpm_go: 1
      mkr_max_zero_p: 0.5
      mkr_do_gsea: true
      mkr_gsea_cut: 0.1
      mkr_gsea_var: z_score
      mkr_custom_genesets:
      - name: mouse_brain
        file: /path/to/file/with/marker/genes.csv
    label_celltypes:
      - labeller: "scprocess"
        model: "human_cns"
        hi_res_cl: "RNA_snn_res.2"
        min_cl_prop: 0.5
        min_cl_size: 100
    train_xgboost:
      annots_f: /path/to/annotations.csv
      ref_tag: my_classifier                      # optional; defaults to xgboost_{full_tag}
      label_map_f: /path/to/label_map.csv
      refine_labels: true
      purity_threshold: 0.65
      n_cells_per_type: 1000
      min_cells_per_type: 20
      min_cells_expressed: 10
      gene_exclude_re: "(lincRNA|lncRNA|pseudogene|antisense)"
      seed: 42
      use_gpu: false
      pass1_subsample: 0.632
      pass1_colsample_bytree: 0.1
      pass1_learning_rate: 0.1
      pass1_nrounds: 300
      pass1_early_stopping: 10
      pass2_colsample_bytree: 0.5
      pass2_learning_rate: 0.05
      pass2_nrounds: 500
      pass2_early_stopping: 10
      gain_threshold: 0.9
      min_genes: 100
      max_genes: 3000
    zoom:
      - /path/to/cell_subset_1_zoom_params.yaml
      - /path/to/cell_subset_2_zoom_params.yaml
      - /path/to/cell_subset_3_zoom_params.yaml
    shiny:
      app_title: Test project shiny app
      email`: testy.mcuser@email.com
      keyword`: cells
      default_gene`: APOE
      n_keep`: 3000
      var_names`: ["Variable 1", "Variable 2", "Variable 3"]
      metadata_combns`:
        - ["Variable 1", "Variable 2"]
      home_md: /path/to/home.md
      annotation_csv`: /path/to/annotation.csv
      cluster_palette: Renoir
      metadata_palettes:
        var1: VanGogh1               # palette name
        var2: ["#FF0000", "#0000FF"]       # explicit colours
        var3:                              # explicit colours with custom ordering
          colours: [red, blue, green]
          values: [young, middle, old]
    resources:
      retries: 3
      n_run_mapping: 8
    ```


#### Required parameters

##### project

* `proj_dir`: absolute path to `workflowr` project directory created with the {{scnew}} function.
* `fastq_dir`: path to directory containing FASTQ files. Should be absolute or relative to `proj_dir`. Exactly one of `fastq_dir` and `arv_uuids` should be specified.
* `arv_uuids`: list of Arvados UUIDs where fastq files are located. Exactly one of `fastq_dir` and `arv_uuids` should be specified.
* `arv_instance`: the name of Arvados instance. Required if `arv_uuids` is defined.
* `full_tag`: full project label, used in output file names.
* `short_tag`: abbreviated project label, used in output directory names.
* `your_name`: author’s name, displayed in HTML outputs.
* `affiliation`: author’s affiliation, displayed in HTML outputs.
* `date_stamp`: start date of the analysis, formatted as `"YYYY-MM-DD"`.
* `sample_metadata`: path to CSV file with sample metadata. Should be absolute or relative to `proj_dir`. Spaces in column names are not allowed. Only required column is `sample_id`; values in `sample_id` should not contain `_R1`/`.R1` and `_R2`/`.R2` strings and should not overlap (a value should not be a subset of any other values).
* `ref_txome`: must match one of the values in the `ref_txome` column of `index_parameters.csv` (created by {{scsetup}}). Required for polyA data; must not be specified for Flex data (use `probe_set` instead).
* `tenx_assay_type`: assay type. Options are `poly_a` and `flex`. When set to `flex`, `probe_set` is required and `ref_txome` must not be specified.
* `probe_set`: probe set to use for 10x Flex data. Required for Flex data. Valid values are `human_v1`, `human_v2`, `mouse_v1`, and `mouse_v2`. Must not be specified for polyA data (use `ref_txome` instead).

#### Optional parameters

##### project

* `tenx_chemistry`: 10x assay configuration. Accepted values are `3v2`, `3v3`, `3v4`, `5v1`, `5v2`, `5v3`, `multiome`, `flexv1` and `flexv2`. `multiome` refers only to gene expression data generated with the 10x multiome kit (ATACseq data is not supported). For Flex data, chemistry specification is redundant as chemistry is derived automatically based on `probe_set`. When not specified for polyA data, chemistry is auto-detected from barcode whitelist overlap.
* `metadata_vars`: A list of column names in the `sample_metadata` file to be used for visualizing the distribution of cell annotations across identified clusters and regions of the low-dimensional embedding.
* `show_arv_uuids`: Whether to display Arvados UUIDs (`arv_uuids`) in the configuration file details box on the index page. If `false`, UUIDs are replaced with "not shown". Defaults to `true`.
* `exclude`: List of all samples that should be excluded from the analysis. Samples can be listed under `pool_id` (if multiplexed) or `sample_id`. 
* `custom_sample_params`: YAML file with optional custom parameters for each pool or sample (custom `tenx_chemistry`, custom `mapping`, custom `ambient` and custom `qc` parameters can be specified for each sample). For example:

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
      cb_posterior_batch_size: 128 # only applicable if cb_version is v.0.3.2
sample_id:
  sample_1:
    qc:
      qc_min_counts: 100
```

##### multiplexing

* `demux_type`: `demux_type` options (default is `none`):
    + `none` if experiment is not multiplexed;
    + `hto` if hto-based demultiplexing of samples should be performed with {{sc}};
    + `ocm` if on-chip multiplexing was used. Requires an `ocm_id` column in sample metadata with values `OB1`-`OB4`.
    + `custom` if demultiplexing results will be used as input to {{sc}}; or
    + `flex` for 10x Flex data where samples are pooled within a library and demultiplexed using probe barcodes. Requires `tenx_assay_type: flex` in the `project` section.
* `fastq_dir`: path to directory containing HTO FASTQ files. Should be absolute or relative to `proj_dir`. If `demux_type` is `hto`, exactly one of `fastq_dir` and `arv_uuids` should be specified.
* `arv_uuids`: list of Arvados UUIDs where HTO fastq files are located. Expects `arv_instance`to be defined. If `demux_type` is `hto`, exactly one of `fastq_dir` and `arv_uuids` should be specified.
* `feature_ref`: path to CSV file with columns `hto_id` and `sequence`. Required if `demux_type` is `hto`.
* `demux_output`: path to CSV file with columns `pool_id`, `sample_id`, `cell_id`. Optional column `class` can be added with values `doublet`, `singlet` or `negative`. Required if `demux_type` is `custom`.
* `seurat_quantile`: equivalent to the `positive.quantile` argument of the `Seurat::HTODemux` function (see [Seurat documentation](https://satijalab.org/seurat/reference/htodemux) for more details). Used if `demux_type` is `hto`.

##### ambient

* `ambient_method`: method for ambient RNA removal. Options are:
    + `decontx_background` (default): runs [`decontX`](https://bioconductor.org/packages/release/bioc/html/celda.html) with an estimated ambient RNA profile as input.
    + `decontx_cluster`: runs `decontX` with a cluster-based background model.
    + `cellbender`: uses [CellBender](https://cellbender.readthedocs.io) for ambient RNA removal.
    + `none`: skips ambient RNA removal;
* `cell_calling`: method for cell calling when `ambient_method` is `none`, `decontx_background`, or `decontx_cluster`. Options are `barcodeRanks` (default) and `emptyDrops`.
* `cb_version`: version of `cellbender` to use if `ambient_method` is set to `cellbender`. Options are `v0.3.2` (default), `v0.3.0'` and `v0.2.0'`.
* `cb_max_prop_kept`: maximum proportion of droplets, relative to `--total-droplets-included`, that `cellbender` can call as cells. Default is `0.9`, meaning samples are excluded if `cellbender` calls more than 90% of `--total-droplets-included` droplets as cells. Applicable only if `ambient_method` is `cellbender`. For more information about the `--total-droplets-included` parameter see [Cellbender documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_learning_rate`: Sets the `--learning-rate` `CellBender` parameter to the specified value; applicable only if `ambient_method` is `cellbender`. Default value is `0.0001`. For more information about this parameter see [Cellbender documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_empty_training_fraction`: Sets the `--empty-drop-training-fraction` `CellBender` parameter to the specified value; applicable only if `ambient_method` is `cellbender`. Default value is `0.2`. Setting this to a lower value (e.g. 0.1 or 0.05) can help if `CellBender` jobs are failing on samples with very few cells. For more information about this parameter see [Cellbender documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_posterior_batch_size`: Value of the `--posterior-batch-size` parameter; applicable only if `ambient_method` is `cellbender` and `cellbender_version` is `v0.3.2`. For more information about this parameter see [Cellbender documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_expected_cells`: forces the `--expected-cells` `Cellbender` parameter to be consistent across all samples; applicable only if `ambient_method` is `cellbender`. For more information about this parameter see [Cellbender documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_total_droplets_included`: forces the `--total-droplets-included` `Cellbender` parameter to be consistent across all samples; applicable only if `ambient_method` is `cellbender`. For more information about this parameter see [Cellbender documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_low_count_threshold`: forces the `--low-count-threshold` `CellBender` parameter to be consistent across all samples; applicable only if `ambient_method` is `cellbender`. For more information about this parameter see [Cellbender documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).


##### qc

* `qc_min_counts`: minimum number of UMIs per cell required to retain the cell.
* `qc_max_counts`: optional maximum number of UMIs allowed per cell; `null` disables the upper limit.
* `qc_min_feats`: minimum number of detected features per cell required to retain the cell.
* `qc_max_feats`: optional maximum number of detected features allowed per cell; `null` disables the upper limit.
* `qc_min_mito`: minimum proportion of mitochondrial reads required to retain the cell.
* `qc_max_mito`: maximum proportion of mitochondrial reads allowed to retain the cell.
* `qc_min_splice`: minimum proportion of spliced reads required to retain the cell.
* `qc_max_splice`: maximum proportion of spliced reads allowed to retain the cell.
* `qc_min_cells`: minimum number of cells required in a sample after QC filtering to retain the sample.
* `dbl_min_feats`: number of features required for each barcode to be included in `scDblFinder` calculations.
* `exclude_mito`: boolean; whether to exclude mitochondrial genes or not.

Configured maxima must be greater than or equal to their corresponding minima. Count and feature maxima are inclusive and are applied during final QC filtering, after `scDblFinder` has run.

##### pb_empties

* `ambient_genes_logfc_thr`: log-fold change (logFC) threshold used to filter the results of the edgeR differential expression test comparing empty droplets to cells.
* `ambient_genes_fdr_thr`: false discovery rate (FDR) threshold used to filter the results of the edgeR differential expression test comparing empty droplets to cells.

##### hvg

* `hvg_method`: options:
    + `sample` - calculate highly variable genes per sample, then calculate combined ranking across samples;
    + `all` - calculate highly variable genes across all cells in the dataset; and
    + `groups` - calculate highly variable genes for each sample group then calculate combined ranking across groups.
* `hvg_metadata_split_var`: if `hvg_method` is `groups`, which variable in `sample_metadata` should be used to define sample groups.
* `hvg_n_hvgs`: number of HVGs to use for PCA
* `hvg_exclude_ambient_genes`: if `true`, genes enriched in empty droplets relative to cells will be excluded from highly variable genes selection.
* `hvg_exclude_from_file`: path to CSV file with genes to be excluded from HVGs. Should be absolute or relative to `proj_dir`. File should contain one column, named either `gene_id` or `symbol`. Values in the column should all be present in reference genome.
* `hvg_chunk_size`: number of genes to use for each chunked matrix.


##### integration

* `int_pca_method`: PCA computation method. `bpcells` (default) streams normalized and scaled expression from disk and is intended for large projects; `scanpy` retains the original in-memory CPU/GPU PCA path. For primary projects, the BPCells path persists preliminary doublet classifications and a clean-cell whitelist before fitting the final clean-cell PCA.
* `int_use_gpu`: whether to use GPU acceleration (`RAPIDS-singlecell`) for integration and clustering steps. Options are `true` (default) or `false`. If GPU is not available, `Scanpy` will be used. For join integration with multiple `int_batch_var` values, setting this to `true` uses CPU Harmony and GPU-accelerated neighbors, Leiden, and UMAP; a warning reports this hybrid execution mode.
* `int_embedding`: which dimensionality reduction method to use for clustering and UMAP, options: `pca` (no batch correction), `harmony` (batch correction). 
* `int_theta`: theta parameter for `Harmony` integration, controlling batch variable mixing.
* `int_batch_var`: variable to use for integration with `Harmony`. Default is `sample_id`; if `demux_type` is set to either `hto` or `custom`, then `pool_id` is an alternative option.
* `int_n_dims`: number of principal components to use for data integration.
* `int_dbl_res`: clustering resolution for identification of additional doublets.
* `int_dbl_cl_prop`: threshold for the proportion of doublets within a cluster. Clusters where the proportion of doublets exceeds this value will be excluded.
* `int_sce_outs`: if `true` H5AD outputs will be converted to `SingleCellExperiment` objects and stored ad RDS files.
* `int_res_ls`: list of resolution values to be used for clustering.
* `int_use_paga`: if `true`, enable Partition-based graph abstraction (PAGA) for trajectory analysis and cell hierarchy inference. A clustering at the specified resolution will be computed for PAGA.
* `int_paga_cl_res`: clustering resolution for PAGA analysis. Must be a value listed in `int_res_ls`. Default is 2 for project, zoom, and join analyses. Only used when `int_use_paga` is `true`. This is independent of `mkr_sel_res`, which defaults to 0.2.

##### cell_cycle

The presence of the optional `cell_cycle` block enables fixed-reference
tricycle projection. If the block is absent, no tricycle files or fields are
created. An empty block (`cell_cycle: {}`) uses the defaults and produces three
cell-level values:

```yaml
cell_cycle: {}
```

Tricycle projections, the pooled origin, diagnostics, and any regression
coefficients are stored in `output/{short_tag}_cell_cycle/`. Per-sample
projection intermediates are stored in its `tricycle/` subdirectory.

* `tricycle_pc1`: first fixed-reference tricycle projection coordinate after
  tricycle mean-centres expression within the biological sample;
* `tricycle_pc2`: second coordinate with the same per-sample centring; and
* `tricycle_theta`: phase angle in radians in `[0, 2*pi)`, recalculated around
  one density-equalized origin estimated from QC-passed cells across the
  project.

The projection coordinates are not integration PCs, and both the per-sample
centring and pooled origin make these values dependent on project composition.
For HTO/custom demultiplexing, doublets without a biological sample assignment
are projected using a clearly recorded run-centred fallback so they can enter
the preliminary doublet pass; they are excluded from origin estimation.

Add the optional `regression` subsection to remove shared cyclic effects from
normalized, scaled HVG expression before both preliminary and final PCA. An
empty subsection (`regression: {}`) uses two harmonics and ridge penalty 0.1.
Regression requires `integration.int_pca_method: bpcells`. `harmonics` is `1`
or `2`, and `ridge_lambda` controls shrinkage of the shared sine/cosine
coefficients. The cyclic columns are centred within sample and scaled globally;
sample-specific slopes are not offered.

```yaml
cell_cycle:
  regression: {}
```

##### marker_genes

* `mkr_sel_res`: selected cluster resolution used for identifying marker genes.
* `mkr_min_cl_size`: minimum number of cells required in a cluster to calculate marker genes for that cluster.
* `mkr_min_cells`: minimum number of cells required in a pseudobulk sample to include it in marker gene calculations.
* `mkr_adjust_project_id`: join-only boolean controlling whether marker-gene
  dispersion and one-versus-rest Treat models include `project_id` as an
  additive blocking factor. The default is `false`, which uses cluster-only
  models and edgeR's faster one-way fitting path. Set it to `true` to account
  for project-level technical differences.
* `mkr_not_ok_re`: regular expression pattern to exclude specific gene types from plots showing marker gene expression.
* `mkr_min_cpm_mkr`: minimum counts per million (CPM) in a cell type required for a gene to be considered a marker gene.
* `mkr_do_gsea`: boolean specifiying whether Gene Set Enrichment Analysis (GSEA) should be performed on marker genes. 
* `mkr_min_cpm_go`: minimum counts per million (CPM) in a cell type required for a gene to be used in GSEA.
* `mkr_max_zero_p`: maximum proportion of pseudobulk samples for a cell type that can have zero counts for a gene to be used in GSEA.
* `mkr_gsea_cut`: False discovery rate (FDR) cutoff for GSEA.
* `mkr_gsea_var`: the statistical measure used for ranking genes for GSEA. Choices are `z_score` (z-score based on signed `log10(FDR)`, the default) or `logFC` (log fold change).
* `mkr_custom_genesets`: a list of custom marker gene sets, each defined by a unique name and associated file path.
    + `name`: a string representing the name of the marker gene set
    + `file`: path to CSV file containing a list of genes in the marker gene set. Must contain column `label` (marker gene category), and `symbol` and/or `ensembl_id`. If not speficied `scprocess` will look for file `$SCPROCESS_DATA_DIR/marker_genes/{name}.csv`

For pseudobulk HVG ranking, cluster-sample profiles containing only one cell are
excluded from the VST matrix; this does not exclude them from marker-gene testing.
The VST uses standard DESeq2 size factors when possible and falls back to the
`poscounts` estimator for sparse matrices in which every gene contains a zero.
Join analyses instead rank variability as the blockwise variance of
TMM-normalised log-CPM across eligible cluster-sample pseudobulks. The join
calculation reads bounded gene blocks from the disk-backed BPCells matrix and
does not construct a dense, whole-matrix VST.

##### label_celltypes

* `labeller`: specifies the method to annotate cell types; options include:
    + `celltypist` uses one of the models available in [`CellTypist`](https://www.celltypist.org/models) for annotation.
    + `scprocess`: use an `XGBoost` classifier for cell type annotation.
* `model`: determines the model to be used based on the selected `labeller`. For list of all available `CellTypist` models see `$SCPROCESS_DATA_DIR/celltypist/celltypist_models.csv`. If `labeller` is set to `scprocess` the value should match the name of a classifier registered with {{scsetup}} (see `$SCPROCESS_DATA_DIR/xgboost/xgboost_models.csv` for a list of available classifiers).
* `hi_res_cl`: name of a column containing high-resolution clustering results. It must follow the pattern `"RNA_snn_res.n"` where `n` should be replaced with one of the values in `int_res_ls`. Default is `"RNA_snn_res.2"`.
* `min_cl_prop`: minimum proportion of cells in a cluster that need to be labeled for that cluster to be labeled.
* `min_cl_size`: minimum number of cells in a cluster required for that cluster to be labeled.
* `label_map`: path to a CSV file mapping fine-grained labels to coarse labels. Must contain columns `fine_label` and `coarse_label`. Should be absolute or relative to `proj_dir`. If not specified, {{sc}} will use the default label map included with the classifier (if available).
* `save_cluster_names_file`: if `true`, generate a `cluster_names_for_shiny_*.csv` file mapping clusters (at the `marker_genes:mkr_sel_res` resolution) to predicted cell type names. This file can be used as the `annotation_csv` in the shiny config. Requires a `marker_genes` block to be configured. Default is `false`.


##### train_xgboost

* `annots_f` (required): path to a CSV file containing cell annotations. Must have columns `cell_id` and `annotation`. Should be absolute or relative to `proj_dir`.
* `ref_tag` (optional): a unique identifier for the classifier.
* `label_map_f`: path to a CSV file mapping fine-grained annotations to coarse labels. Must have columns `annotation` and `coarse_label`. If provided, the trained classifier will also be evaluated on coarse labels and a label map file will be included in the classifier outputs. Should be absolute or relative to `proj_dir`.
* `refine_labels`: if `true` (default), labels are refined by cluster majority voting — clusters where the majority of cells share the same annotation are assigned that label, which helps correct noisy annotations.
* `purity_threshold`: minimum proportion of cells in a cluster that must share the same annotation for the cluster to be assigned that label during label refinement. Default is `0.65`.
* `n_cells_per_type`: target number of cells per cell type after downsampling to balance the training set. Default is `1000`.
* `min_cells_per_type`: minimum number of cells required to retain a cell type in the training set. Cell types with fewer cells are excluded. Default is `20`.
* `min_cells_expressed`: minimum number of cells in which a gene must be detected (non-zero counts) for the gene to be included as a feature. Default is `10`.
* `gene_exclude_re`: regular expression pattern to exclude specific gene types from training. Default is `"(lincRNA|lncRNA|pseudogene|antisense)"`; this default allows the classifier to be used on probe-based Flex data as well as reverse-transcription-based data.
* `seed`: random seed for reproducibility. Default is `42`.
* `use_gpu`: whether to use GPU acceleration for `XGBoost` training. Default is `false`.
* `gain_threshold`: cumulative gain fraction at which to stop selecting genes. Default is `0.9` (top 90% of total gain).
* `min_genes`: minimum number of genes to select, regardless of gain threshold. Default is `100`.
* `max_genes`: maximum number of genes to select. Default is `3000`.

??? note "XGBoost training parameters"

    The classifier is trained in two passes. The first pass trains on all genes to compute feature importance scores. Genes are then selected based on cumulative gain, and the second pass trains the final model on the selected gene set.

    **Pass 1 parameters** (broad feature exploration):

    * `pass1_subsample`: fraction of training rows sampled per boosting round. Default is `0.632`.
    * `pass1_colsample_bytree`: fraction of features (genes) sampled per tree. Default is `0.1`.
    * `pass1_learning_rate`: learning rate (step size shrinkage). Default is `0.1`.
    * `pass1_nrounds`: maximum number of boosting rounds. Default is `300`.
    * `pass1_early_stopping`: stop training if validation performance does not improve for this many rounds. Default is `10`.

    **Pass 2 parameters** (refined training on selected genes):

    * `pass2_colsample_bytree`: fraction of selected genes sampled per tree. Default is `0.5`.
    * `pass2_learning_rate`: learning rate. Default is `0.05`.
    * `pass2_nrounds`: maximum number of boosting rounds. Default is `500`.
    * `pass2_early_stopping`: stop training if validation performance does not improve for this many rounds. Default is `10`.

##### shiny

* `app_tag`: machine-safe identifier used for the deployment directory, generated data-file prefixes, sentinel, and build log. It may contain letters, numbers, `.`, `_`, and `-`. Defaults to `project.short_tag` for a standard project and `join.name` for a join. Main apps are written to `public/shiny_{app_tag}`.
* `app_title`: title displayed in the Shiny app header. Defaults to `short_tag`.
* `email`: contact email shown in the app footer.
* `keyword`: short word used in plot axis labels and descriptions (e.g. `"cells"`, `"nuclei"`). Default is `"cells"`.
* `default_gene`: gene symbol displayed by default in the "Explore Genes" tab.
* `n_keep`: number of cells retained in the subsampled UMAP shown in the app. Default is 30000.
* `metadata_vars`: metadata variables to show in the app. If not specified, uses the project-level `metadata_vars`. Use this to show a different set of variables in the Shiny app than in the main pipeline.
* `metadata_labels`: display labels for metadata variables. Keys are column names, values are display labels. Variables not listed keep their column name as the label. Example:
    ```yaml
    shiny:
      app_tag: test_project
      metadata_labels:
        brainregion: "brain region"
        condition: "treatment group"
    ```
* `var_names`: _(deprecated, use `metadata_labels` instead)_ display names for `metadata_vars` columns (same order). Defaults to `metadata_vars` values.
* `metadata_combns`: list of metadata variable pairs to display as combined groupings. Each element should be a two-element list of variable names from `metadata_vars` (use column names, not display labels).
* `home_md`: path to a Markdown file used as the landing page content. Should be absolute or relative to `proj_dir`.
* `annotation_csv`: path to a CSV file with columns `cluster`, `cluster_name`, and optionally `colour`, defining display names, order, and colours for clusters. Absolute or relative to `proj_dir`. During app configuration, its cluster values must exactly match the full set of clusters in the built Shiny data, including rare clusters without marker-gene results. Missing or extra clusters, duplicated cluster IDs or display names, and missing or blank cluster IDs or names cause configuration to stop with an error; this prevents annotations from a different clustering resolution or analysis run being applied silently. UMAP subsampling retains at least one cell from every cluster.
* `cluster_palette`: name of colour palette applied to clusters when `annotation_csv` is not provided. Accepts any name from the Supported colour palettes list.
* `metadata_palettes`: per-variable colour palette configuration. Each key is a metadata variable name; the value can be:
    - a palette name (string)
    - an explicit colour list (array)
    - an object with optional `palette`, `colours`, and `values` keys (the `values` list controls level ordering; defaults to frequency order)

Shiny deployment has two stages. `build_shiny_data` (or
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
Shiny data. The first run after upgrading an existing deployment to this split
creates the new data-stage sentinel and may rebuild the data once.

The cluster UMAPs include a **Repel cluster labels** option. Repelled label
positions are computed from the full UMAP during the app build: labels are moved
into low-density space, given clearance from cells, and allocated globally so
their numbered circles do not overlap one another. Rebuild the Shiny target to
refresh these positions after changing the integration or clustering.

The adjacent cluster legend sizes itself from the longest annotation while
favouring space for the UMAP panel; labels may extend into a small unclipped
margin rather than reserving a large empty legend area.

The sample metadata CSV is an explicit input to the Shiny data stage. Editing
descriptive values therefore marks the generated Shiny data for rebuilding. For a join, first
run the normal join target so its joint sample metadata snapshot is refreshed, then
run `scprocess join CONFIG -r shiny` to rebuild the app.


??? info "Supported colour palettes"

    The following palette names are accepted anywhere a `cluster_palette` or `metadata_palettes` entry expects a name:

    | Source | Names |
    |--------|-------|
    | **scprocess built-in** | `nice_cols` (42 colours) |
    | **MetBrewer** | `Archambault`, `Austria`, `Benedictus`, `Cassatt1`, `Cassatt2`, `Cross`, `Degas`, `Demuth`, `Derain`, `Egypt`, `Gauguin`, `Greek`, `Hiroshige`, `Hokusai1`, `Hokusai2`, `Hokusai3`, `Homer1`, `Homer2`, `Ingres`, `Isfahan1`, `Isfahan2`, `Java`, `Johnson`, `Juarez`, `Kandinsky`, `Klimt`, `Lakota`, `Manet`, `Monet`, `Moreau`, `Morgenstern`, `Nattier`, `Navajo`, `NewKingdom`, `Nizami`, `OKeeffe1`, `OKeeffe2`, `Paquin`, `Peru1`, `Peru2`, `Pillement`, `Pissaro`, `Redon`, `Renoir`, `Signac`, `Tam`, `Tara`, `Thomas`, `Tiepolo`, `Troy`, `Tsimshian`, `VanGogh1`, `VanGogh2`, `VanGogh3`, `Veronese`, `Wissing` |
    | **RColorBrewer** | Qualitative: `Accent`, `Dark2`, `Paired`, `Pastel1`, `Pastel2`, `Set1`, `Set2`, `Set3`. Sequential/diverging: `Blues`, `BuGn`, `BuPu`, `GnBu`, `Greens`, `Greys`, `Oranges`, `OrRd`, `PuBu`, `PuBuGn`, `PuRd`, `Purples`, `RdPu`, `Reds`, `YlGn`, `YlGnBu`, `YlOrBr`, `YlOrRd`, `BrBG`, `PiYG`, `PRGn`, `PuOr`, `RdBu`, `RdGy`, `RdYlBu`, `RdYlGn`, `Spectral`. |
    | **ggsci** | `npg`, `aaas`, `nejm`, `lancet`, `jama`, `jco`, `ucscgb`, `d3` / `d3_20`, `d3_10`, `d3_20b`, `d3_20c`, `igv`, `cosmic`, `simpsons`, `rickandmorty`, `futurama`, `tron`, `startrek`, `uchicago`, `frontiers`, `flatui`, `bootstrap` |


##### zoom

In this section, users can provide multiple YAML files, each specifying parameters for repeating certain steps of {{sc}} on a subset of cells. Some parameters in the YAML file inherit their definitions from the primary {{sc}} configuration file, including `qc_min_cells`, `qc_min_counts`, `qc_max_counts`, `qc_min_feats`, `qc_max_feats`, `qc_min_mito`, `qc_max_mito`, `qc_min_splice`, `qc_max_splice`, `hvg_method`, `hvg_metadata_split_var`, `hvg_n_hvgs`, `hvg_chunk_size`, `hvg_exclude_ambient_genes`, `hvg_exclude_from_file`, `ambient_genes_logfc_thr`, `ambient_genes_fdr_thr`, `int_pca_method`, `int_use_gpu`, `int_embedding`, `int_n_dims`, `int_theta`, `int_res_ls`, `int_use_paga`, `int_paga_cl_res`, `mkr_sel_res`, `mkr_min_cl_size`, `mkr_min_cells`, `mkr_not_ok_re`, `mkr_min_cpm_mkr`, `mkr_min_cpm_go`, `mkr_max_zero_p`, `mkr_do_gsea`, `mkr_gsea_cut`, `mkr_gsea_var`,`mkr_custom_genesets`, all [shiny app parameters](#shiny) (`app_title` defaults to `name`), and all [XGBoost training parameters](#train_xgboost).

Additional parameters include:

* `name`: name of cell subset to be analysed. 
* `labels_source`: specifies how a cell subset is defined (required). Options include:
    - `scprocess`: labels assigned by the `XGBoost` classifier (using rule `label_celltypes`)
    - `celltypist`: labels assigned by `CellTypist`(using the rule `label_celltypes`)
    - `clusters`: labels based on clustering results obtained with {{sc}}
    - `custom`: user-defined cell type annotations
* `model`: required if `labels_source` is `scprocess` or `celltypist`. 
* `sel_labels`: a list of all labels that define cell types/clusters to be included in subclustering (required).
+ `labels_col`: name of column that contains cell type/cluster labels.
* `save_subset_sces`: whether to create `SingleCellExperiment` objects containing cells that have been assigned one of the values in `sel_labels`; default is `false`.
* `save_subset_anndata`: whether to create H5AD files containing cells that have been assigned one of the values in `sel_labels`; defaults is `true`.
* `save_cluster_names_file`: whether to generate a resolution-specific `cluster_names_for_shiny_*.csv` for the zoom clusters. Names are assigned by majority vote from the retained cells' `predicted_label_naive` values and made unique when the same predicted subtype is assigned to multiple clusters. This requires a `celltypist` or `scprocess` labels source and can be supplied as `shiny.annotation_csv`. Default is `false`.
* `custom_labels_f`: required if `labels_source` is set to `custom`; path to CSV file with columns `sample_id`, `cell_id` and `label`.
* `exclude`: optionally excludes batches from this zoom only. Use `sample_id` or `pool_id` to match the project's batch variable, for example `exclude: {sample_id: [sample_a, sample_b]}`. This does not exclude the batches from the standard analysis or other zooms.

If a value in `sel_labels` is absent from the selected labels column, {{sc}} reports a warning and continues with the labels that are present. At least one selected label must be present, and the retained subset must contain enough cells for downstream QC, HVG detection, and integration.

Zoom configs also support optional cell-level QC thresholds to apply stricter filtering to cells that already passed the main pipeline's QC. When set, cells in the zoom subset that do not meet these thresholds are removed before HVG detection and integration. All thresholds default to `null` (no additional filtering):

The zoom report always shows QC distributions. Without zoom-specific thresholds it contains one distributions view; when thresholds are configured it also shows the before/after filtering views.

When a zoom is defined from `celltypist` or `scprocess` labels, the zoom report also reuses the retained cells' naive `label_celltypes` predictions. It re-aggregates those predictions by majority vote against the zoom analysis high-resolution clusters, using the matching parent `label_celltypes` entry's `hi_res_cl` (default `RNA_snn_res.2`) and `min_cl_prop`. The report shows naive predicted labels against those high-resolution zoom clusters, zoom-aggregated labels over the zoom UMAP, and naive-versus-aggregated label totals. Coarse-label panels are included when the source label output contains coarse predictions. The zoom workflow does not rerun the classifier.

* `qc_min_counts`: minimum UMI counts per cell. Only effective if stricter than the main pipeline threshold. Default: not set.
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

This section allows users to adjust the resource requirements for specific `Snakemake rule`s. This is especially useful when a step/rule fails on a cluster due to insufficient memory or runtime limits. By specifying the parameters below, users can fine-tune these settings for their pipeline:

* `gb_[rule_name]`: specifies the maximum memory (in GB) requested for running a specific rule. `rule_name` should be replaced with an {{sc}} rule name. This value applies for the entire job, not per thread.
* `mins_[rule_name]`: specifies the maximum runtime (in minutes) requested for running a specific rule. `rule_name` should be replace with an {{sc}} rule name.

By default, {{sc}} predicts resource requirements using linear models fitted on benchmark data from ~96 projects. The prediction formula is `max(floor, intercept + slope × x)`, where `x` is a rule-specific predictor such as input file size or number of batches. Five per-sample rules have separate model parameters for 10x v3 vs v4 chemistry; for auto-detected chemistry, the resolved chemistry is read from the mapping output.

User-specified values always override the model predictions. If a rule has no fitted model and no user override, it falls back to a schema default.

Additional parameters include:

* `retries`: number of times to retry running a specific rule in {{sc}} if it fails. For each attempt the initial memory requested for the rule is multiplied by `1.5**(attempt - 1)`. Useful for when {{sc}} is ran on a [cluster](setup.md#cluster-setup).
* `n_run_mapping`: number of threads requested for running the mapping step if `tenx_assay_type` is `poly_a`. Default is 8.
* `n_run_mapping_flex`: number of threads requested for running the mapping step if `tenx_assay_type` is `flex`. Default is 8.

??? note "Detailed information about resource parameters"

    * `gb_build_hto_index`: maximum memory required (in GB) for rule `build_hto_index`.
    * `gb_run_mapping`: maximum memory required (in GB) for rule `run_mapping`.
    * `gb_run_mapping_flex`: maximum memory required (in GB) for rule `run_mapping_flex`. Default is 12 GB.
    * `gb_save_alevin_flex_to_h5`: maximum memory required (in GB) for rule `save_alevin_flex_to_h5`. Default is 32 GB.
    * `gb_run_mapping_hto`: maximum memory required (in GB) for rule `run_mapping_hto`.
    * `gb_save_alevin_to_h5`: maximum memory required (in GB) for rule `save_alevin_to_h5`.
    * `gb_make_hto_sce_objects`: maximum memory required (in GB) for rule `make_hto_sce_objects`.
    * `gb_save_alevin_hto_to_h5`: maximum memory required (in GB) for rule `save_alevin_hto_to_h5`.
    * `gb_run_cellbender`: maximum memory required (in GB) for rule `run_cellbender`.
    * `gb_run_decontx`: maximum memory required (in GB) for rule `run_decontx`.
    * `gb_run_cell_celling`: maximum memory required (in GB) for rule `run_cell_celling`.
    * `gb_get_barcode_qc_metrics`: maximum memory required (in GB) for rule `get_barcode_qc_metrics`.
    * `gb_run_qc_one_run`: maximum memory required (in GB) for rule `run_qc_one_run`.
    * `gb_merge_qc`: maximum memory required (in GB) for rule `merge_qc`.
    * `gb_merge_rowdata`: maximum memory required (in GB) for rule `merge_rowdata`.
    * `gb_get_qc_sample_statistics`: maximum memory required (in GB) for rule `get_qc_sample_statistics`.
    * `gb_make_one_pb_empty`: maximum memory required (in GB) for rule `make_one_pb_empty`.
    * `gb_merge_pb_empty`: maximum memory required (in GB) for rule `merge_pb_empty`.
    * `gb_make_one_pb_cells`: maximum memory required (in GB) for rule `make_one_pb_cells`.
    * `gb_merge_pb_cells`: maximum memory required (in GB) for rule `merge_pb_cells`.
    * `gb_calculate_ambient_genes`: maximum memory required (in GB) for rule `calculate_ambient_genes`.
    * `gb_make_tmp_csr_matrix`: maximum memory required (in GB) for rule `make_tmp_csr_matrix`.
    * `gb_get_stats_for_std_variance_for_sample`: maximum memory required (in GB) for rule `get_stats_for_std_variance_for_sample`.
    * `gb_merge_sample_std_var_stats`: maximum memory required (in GB) for rule `merge_sample_std_var_stats`.
    * `gb_get_mean_var_for_group`: maximum memory required (in GB) for rule `get_mean_var_for_group`.
    * `gb_get_estimated_variances`: maximum memory required (in GB) for rule `get_estimated_variances`.
    * `gb_get_stats_for_std_variance_for_group`: maximum memory required (in GB) for rule `get_stats_for_std_variance_for_group`.
    * `gb_get_highly_variable_genes`: maximum memory required (in GB) for rule `get_highly_variable_genes`.
    * `gb_create_hvg_matrix`: maximum memory required (in GB) for rule `create_hvg_matrix`.
    * `gb_create_doublets_hvg_matrix`: maximum memory required (in GB) for rule `create_doublets_hvg_matrix`.
    * `gb_run_integration`: maximum memory required (in GB) for rule `run_integration`.
    * `gb_make_clean_h5ads`: maximum memory required (in GB) for rule `make_clean_h5ads`.
    * `gb_convert_h5ad_to_sce`: maximum memory required (in GB) for rule `convert_h5ad_to_sce`.
    * `gb_run_marker_genes`: maximum memory required (in GB) for rule `run_marker_genes`.
    * `gb_run_fgsea`: maximum memory required (in GB) for rule `run_fgsea`.
    * `gb_render_html_mapping`: maximum memory required (in GB) for rule `render_html_mapping`.
    * `gb_render_html_multiplexing`: maximum memory required (in GB) for rule `render_html_multiplexing`.
    * `gb_render_html_ambient`: maximum memory required (in GB) for rule `render_html_ambient`.
    * `gb_render_html_qc`: maximum memory required (in GB) for rule `render_html_qc`.
    * `gb_render_html_hvgs`: maximum memory required (in GB) for rule `render_html_hvgs`.
    * `gb_render_html_integration`: maximum memory required (in GB) for rule `render_html_integration`.
    * `gb_render_html_marker_genes`: maximum memory required (in GB) for rule `render_html_marker_genes`.
    * `gb_render_html_label_celltypes`: maximum memory required (in GB) for rule `render_html_label_celltypes`.
    * `gb_run_celltypist`: maximum memory required (in GB) for rule `run_celltypist`.
    * `gb_run_scprocess_labeller`: maximum memory required (in GB) for rule `run_scprocess_labeller`.
    * `gb_merge_labels`: maximum memory required (in GB) for rule `merge_labels`.
    * `gb_train_xgboost_train`: maximum memory required (in GB) for rule `train_xgboost_train`.
    * `gb_train_xgboost_predict`: maximum memory required (in GB) for rule `train_xgboost_predict`.
    * `gb_train_xgboost_aggregate`: maximum memory required (in GB) for rule `train_xgboost_aggregate`.
    * `gb_render_html_train_xgboost`: maximum memory required (in GB) for rule `render_html_train_xgboost`.
    * `gb_zoom_get_sample_statistics`: maximum memory required (in GB) for rule `zoom_get_sample_statistics`.
    * `gb_zoom_make_one_pb_cells`: maximum memory required (in GB) for rule `zoom_make_one_pb_cells`.
    * `gb_zoom_calculate_ambient_genes`: maximum memory required (in GB) for rule `zoom_calculate_ambient_genes`.
    * `gb_zoom_make_hvg_df`: maximum memory required (in GB) for rule `zoom_make_hvg_df`.
    * `gb_zoom_make_tmp_csr_matrix`: maximum memory required (in GB) for rule `zoom_make_tmp_csr_matrix`.
    * `gb_zoom_get_stats_for_std_variance_for_sample`: maximum memory required (in GB) for rule `zoom_get_stats_for_std_variance_for_sample`.
    * `gb_zoom_get_mean_var_for_group`: maximum memory required (in GB) for rule `zoom_get_mean_var_for_group`.
    * `gb_zoom_merge_group_mean_var`: maximum memory required (in GB) for rule `zoom_merge_group_mean_var`.
    * `gb_zoom_get_estimated_variances`: maximum memory required (in GB) for rule `zoom_get_estimated_variances`.
    * `gb_zoom_get_stats_for_std_variance_for_group`: maximum memory required (in GB) for rule `zoom_get_stats_for_std_variance_for_group`.
    * `gb_zoom_merge_stats_for_std_variance`: maximum memory required (in GB) for rule `zoom_merge_stats_for_std_variance`.
    * `gb_zoom_get_highly_variable_genes`: maximum memory required (in GB) for rule `zoom_get_highly_variable_genes`.
    * `gb_zoom_create_hvg_matrix`: maximum memory required (in GB) for rule `zoom_create_hvg_matrix`.
    * `gb_zoom_run_integration`: maximum memory required (in GB) for rule `zoom_run_integration`.
    * `gb_zoom_run_marker_genes`: maximum memory required (in GB) for rule `zoom_run_marker_genes`.
    * `gb_zoom_run_fgsea`: maximum memory required (in GB) for rule `zoom_run_fgsea`.
    * `gb_zoom_make_subsets`: maximum memory required (in GB) for rule `zoom_make_subsets`.
    * `gb_render_html_zoom`: maximum memory required (in GB) for rule `render_html_zoom`.
    * `mins_build_hto_index`: maximum runtime required (in minutes) for rule `build_hto_index`.
    * `mins_run_mapping`: maximum runtime required (in minutes) for rule `run_mapping`.
    * `mins_run_mapping_flex`: maximum runtime required (in minutes) for rule `run_mapping_flex`. Default is 180 minutes.
    * `mins_save_alevin_flex_to_h5`: maximum runtime required (in minutes) for rule `save_alevin_flex_to_h5`. Default is 10 minutes.
    * `mins_run_mapping_hto`: maximum runtime required (in minutes) for rule `run_mapping_hto`.
    * `mins_save_alevin_to_h5`: maximum runtime required (in minutes) for rule `save_alevin_to_h5`.
    * `mins_make_hto_sce_objects`: maximum runtime required (in minutes) for rule `make_hto_sce_objects`.
    * `mins_save_alevin_hto_to_h5`: maximum runtime required (in minutes) for rule `save_alevin_hto_to_h5`.
    * `mins_run_cellbender`: maximum runtime required (in minutes) for rule `run_cellbender`.
    * `mins_run_decontx`: maximum runtime required (in minutes) for rule `run_decontx`.
    * `mins_run_cell_celling`: maximum runtime required (in minutes) for rule `run_cell_celling`.
    * `mins_get_barcode_qc_metrics`: maximum runtime required (in minutes) for rule `get_barcode_qc_metrics`.
    * `mins_run_qc_one_run`: maximum runtime required (in minutes) for rule `run_qc_one_run`.
    * `mins_merge_qc`: maximum runtime required (in minutes) for rule `merge_qc`.
    * `mins_merge_rowdata`: maximum runtime required (in minutes) for rule `merge_rowdata`.
    * `mins_get_qc_sample_statistics`: maximum runtime required (in minutes) for rule `get_qc_sample_statistics`.
    * `mins_make_one_pb_empty`: maximum runtime required (in minutes) for rule `make_one_pb_empty`.
    * `mins_merge_pb_empty`: maximum runtime required (in minutes) for rule `merge_pb_empty`.
    * `mins_make_one_pb_cells`: maximum runtime required (in minutes) for rule `make_one_pb_cells`.
    * `mins_merge_pb_cells`: maximum runtime required (in minutes) for rule `merge_pb_cells`.
    * `mins_calculate_ambient_genes`: maximum runtime required (in minutes) for rule `calculate_ambient_genes`.
    * `mins_make_tmp_csr_matrix`: maximum runtime required (in minutes) for rule `make_tmp_csr_matrix`.
    * `mins_get_stats_for_std_variance_for_sample`: maximum runtime required (in minutes) for rule `get_stats_for_std_variance_for_sample`.
    * `mins_merge_sample_std_var_stats`: maximum runtime required (in minutes) for rule `merge_sample_std_var_stats`.
    * `mins_get_mean_var_for_group`: maximum runtime required (in minutes) for rule `get_mean_var_for_group`.
    * `mins_get_estimated_variances`: maximum runtime required (in minutes) for rule `get_estimated_variances`.
    * `mins_get_stats_for_std_variance_for_group`: maximum runtime required (in minutes) for rule `get_stats_for_std_variance_for_group`.
    * `mins_get_highly_variable_genes`: maximum runtime required (in minutes) for rule `get_highly_variable_genes`.
    * `mins_create_hvg_matrix`: maximum runtime required (in minutes) for rule `create_hvg_matrix`.
    * `mins_create_doublets_hvg_matrix`: maximum runtime required (in minutes) for rule `create_doublets_hvg_matrix`.
    * `mins_run_integration`: maximum runtime required (in minutes) for rule `run_integration`.
    * `mins_make_clean_h5ads`: maximum runtime required (in minutes) for rule `make_clean_h5ads`.
    * `mins_convert_h5ad_to_sce`: maximum runtime required (in minutes) for rule `convert_h5ad_to_sce`.
    * `mins_run_marker_genes`: maximum runtime required (in minutes) for rule `run_marker_genes`.
    * `mins_run_fgsea`: maximum runtime required (in minutes) for rule `run_fgsea`.
    * `mins_render_html_mapping`: maximum runtime required (in minutes) for rule `render_html_mapping`.
    * `mins_render_html_multiplexing`: maximum runtime required (in minutes) for rule `render_html_multiplexing`.
    * `mins_render_html_ambient`: maximum runtime required (in minutes) for rule `render_html_ambient`.
    * `mins_render_html_qc`: maximum runtime required (in minutes) for rule `render_html_qc`.
    * `mins_render_html_hvgs`: maximum runtime required (in minutes) for rule `render_html_hvgs`.
    * `mins_render_html_integration`: maximum runtime required (in minutes) for rule `render_html_integration`.
    * `mins_render_html_marker_genes`: maximum runtime required (in minutes) for rule `render_html_marker_genes`.
    * `mins_render_html_label_celltypes`: maximum runtime required (in minutes) for rule `render_html_label_celltypes`.
    * `mins_run_celltypist`: maximum runtime required (in minutes) for rule `run_celltypist`.
    * `mins_run_scprocess_labeller`: maximum runtime required (in minutes) for rule `run_scprocess_labeller`.
    * `mins_merge_labels`: maximum runtime required (in minutes) for rule `merge_labels`.
    * `mins_train_xgboost_train`: maximum runtime required (in minutes) for rule `train_xgboost_train`.
    * `mins_train_xgboost_predict`: maximum runtime required (in minutes) for rule `train_xgboost_predict`.
    * `mins_train_xgboost_aggregate`: maximum runtime required (in minutes) for rule `train_xgboost_aggregate`.
    * `mins_render_html_train_xgboost`: maximum runtime required (in minutes) for rule `render_html_train_xgboost`.
    * `mins_zoom_get_sample_statistics`: maximum runtime required (in minutes) for rule `zoom_get_sample_statistics`.
    * `mins_zoom_make_one_pb_cells`: maximum runtime required (in minutes) for rule `zoom_make_one_pb_cells`.
    * `mins_zoom_calculate_ambient_genes`: maximum runtime required (in minutes) for rule `zoom_calculate_ambient_genes`.
    * `mins_zoom_make_hvg_df`: maximum runtime required (in minutes) for rule `zoom_make_hvg_df`.
    * `mins_zoom_make_tmp_csr_matrix`: maximum runtime required (in minutes) for rule `zoom_make_tmp_csr_matrix`.
    * `mins_zoom_get_stats_for_std_variance_for_sample`: maximum runtime required (in minutes) for rule `zoom_get_stats_for_std_variance_for_sample`.
    * `mins_zoom_get_mean_var_for_group`: maximum runtime required (in minutes) for rule `zoom_get_mean_var_for_group`.
    * `mins_zoom_merge_group_mean_var`: maximum runtime required (in minutes) for rule `zoom_merge_group_mean_var`.
    * `mins_zoom_get_estimated_variances`: maximum runtime required (in minutes) for rule `zoom_get_estimated_variances`.
    * `mins_zoom_get_stats_for_std_variance_for_group`: maximum runtime required (in minutes) for rule `zoom_get_stats_for_std_variance_for_group`.
    * `mins_zoom_merge_stats_for_std_variance`: maximum runtime required (in minutes) for rule `zoom_merge_stats_for_std_variance`.
    * `mins_zoom_get_highly_variable_genes`: maximum runtime required (in minutes) for rule `zoom_get_highly_variable_genes`.
    * `mins_zoom_create_hvg_matrix`: maximum runtime required (in minutes) for rule `zoom_create_hvg_matrix`.
    * `mins_zoom_run_integration`: maximum runtime required (in minutes) for rule `zoom_run_integration`.
    * `mins_zoom_run_marker_genes`: maximum runtime required (in minutes) for rule `zoom_run_marker_genes`.
    * `mins_zoom_run_fgsea`: maximum runtime required (in minutes) for rule `zoom_run_fgsea`.
    * `mins_zoom_make_subsets`: maximum runtime required (in minutes) for rule `zoom_make_subsets`.
    * `mins_render_html_zoom`: maximum runtime required (in minutes) for rule `render_html_zoom`.
    




## {{scjoin}} { #scprocess-join }

**Description**: Integrate multiple completed {{sc}} projects.

**Parameters**:

* `configfile` (positional): path to a configuration YAML file.
* `-r`/`--rule` (optional): `all` (default), `shiny` (build the Shiny app from join outputs) or `train_xgboost`(train an XGBoost classifier from join outputs)
* `-n`/`--dry-run` (optional): perform a trial run which lists all steps that {{scjoin}} would do and does not create any new files. Helpful for checking input files and parameters.
* `-E`/`--extraagrs` (optional): list of additional arguments to pass to `Snakemake`. Refer to [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli).
* `--unlock` (optional): unlock the directory if a previous run was interrupted.

### configuration file

The join configuration reuses the `hvg`, `integration`, `marker_genes`, `label_celltypes`, `train_xgboost`, `shiny`, and `resources` sections from a project configuration where supported by the join schema. Use the current generated template and schema for exact fields and defaults. The join must specify a `ref_txome`, or a `probe_set` for Flex data, that matches every source project; only compatible projects can be joined.

Join integration may use arrays for `int_batch_var` and `int_theta` to configure
multiple Harmony covariates. When both values are arrays, they must have the same
length so that each batch variable has a corresponding theta value.

Additional parameters include:
* `name`: short name; output directories are {name}_join and {name}_marker_genes
* `metadata_vars` (optional): list of column names from the source projects' `sample_metadata` CSV files to carry through into the joint sample metadata and presentation outputs. Each variable must exist in at least one project's metadata file, but does not need to be present in all projects. Missing values are filled with `NA`.
* `projects`: mapping containing at least two source projects. Each entry requires `config`, the path to a completed project config, and may specify `zoom_name` to join a completed zoom instead of the full project.
* `int_pca_method`: PCA computation method. Options: `bpcells` (default) uses disk-backed SVD via BPCells (R), suitable for very large datasets (>1M cells) without GPU memory limits; `scanpy` uses the standard in-memory PCA on GPU/CPU (original behaviour).
* Join marker-gene testing stores the pseudobulk output as a BPCells directory,
  keeps the combined counts disk-backed, and uses the streamed `edger.bp` edgeR
  Treat backend. The work is split into `join_make_pseudobulks`,
  `join_prepare_pseudobulks`, `join_calc_hvgs`, and `join_marker_genes`.
  By default, marker-gene models use cluster alone. Set
  `marker_genes.mkr_adjust_project_id: true` to include `project_id` as an
  additive blocking factor in both dispersion estimation and one-versus-rest
  Treat tests. Project adjustment accounts for project-level technical
  differences but uses edgeR's slower general-GLM fitting path.
  Preparation writes small column and gene manifests containing the retained
  pseudobulks, retained genes, and TMM library-size metadata; it does not
  duplicate or mutate the BPCells count store. Once preparation finishes, the
  blockwise, disk-backed variability ranking and marker testing can run in
  parallel. Their default resource requests are 8–16 GB and 30–90 minutes,
  and can be overridden with the corresponding `gb_join_*` and
  `mins_join_*` resource keys.
  `resources.n_join_marker_genes` controls the threads requested by
  `join_marker_genes` and the number of streamed edgeR workers; its default is
  `8`. BLAS and OpenMP are restricted to one thread inside each worker to avoid
  nested parallelism.
  The rule environment installs the audited BPCells and `edger.bp` R 4.5
  builds from the public
  [edger.bp GitHub Conda channel](https://github.com/wmacnair/edger.bp/tree/conda-channel);
  creating this environment for the first time requires access to
  `raw.githubusercontent.com`.
  Standard project and zoom marker-gene runs retain their existing RDS,
  dense-edgeR, and DESeq2 VST paths.

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
  # probe_set: human_v1           # use instead of ref_txome for flex projects
  # tenx_assay_type: flex         # set to 'flex' when joining flex projects
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
    save_cluster_names_file: true
train_xgboost:
  annots_f: /path/to/annotations.csv.gz
  ref_tag: my_classifier
shiny:
  app_tag: my_join
  app_title: My Joint Analysis
```

**Notes on `label_celltypes` in join:**

* When `label_celltypes` is configured, it runs as part of the default `scprocess join` target.
* If a source project already has label_celltypes outputs for the same `labeller`/`model`, naive predictions are reused instead of re-running the labeller — only projects without existing labels trigger fresh runs.
* `save_cluster_names_file: true` generates a `cluster_names_for_shiny_*.csv` at the `marker_genes:mkr_sel_res` resolution. This file can be used as `annotation_csv` in the `shiny:` section.
