# Reference

## {{scsetup}} { #scprocess-setup }

**Description**: Download all data required for {{sc}} and index reference transcriptomes for `simpleaf`.

**Parameters**:

* `-n`/`--dry-run` (optional): list the setup steps without executing them.
* `--create-envs` (optional): only create the Conda environments used by {{scrun}} and {{scjoin}}. Does not need `--rangerurl`. With `-n`, list the environments without creating them.
* `-E`/`--extraargs` (optional): list of additional arguments to pass to `Snakemake`. Refer to [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for a detailed explanation of available command-line options.
* `-c`/`--rangerurl`:  download link for Cell Ranger (v9.0.0 or higher) available on the [10x Genomics CellRanger download & installation page](https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions); required when Cell Ranger resources have not yet been installed.

### configuration file

The command requires a configuration file named `scprocess_setup.yaml` located in {{sc}} data directory (for instructions on how to set up the {{sc}} data directory see the [Getting started](setup.md#scprocess-data-directory-setup) section). In this file, the user can specify parameters that are used across all {{sc}} projects, such as HPC configuration and reference genomes that will be made available for {{sc}}. For example:

```yaml
user:
  profile_template: slurm_default
  profile_name: slurm_my_cluster
  your_name: Testy McUser
  affiliation: where you work
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

* `profile_template`: the name of the bundled HPC profile template to copy from the {{sc}} _profiles_ directory into `$SCPROCESS_DATA_DIR/profiles`. Bundled templates include `slurm_default`, `lsf_default`, `pbs_default`, `sge_default`, and `htcondor_default`. Exactly one of `profile_template` and `local_cores` should be specified.
* `profile_name` (optional): name of the active profile directory under `$SCPROCESS_DATA_DIR/profiles`. Defaults to the value of `profile_template`.
* `local_cores`: number of CPU cores available for local execution (see [Snakemake documentation](https://snakemake.readthedocs.io/en/v9.8.0/executing/cli.html) for more details). Exactly one of `profile_template` and `local_cores` should be specified.
* `your_name` (optional): author's name. If specified it will be used in the configuration file for new projects created with the `scprocess newproj -c` command.
* `affiliation` (optional): author's affiliation. If specified it will be used in the configuration file for new projects created with the `scprocess newproj -c` command.
* `int_use_gpu` (optional): whether to use GPU acceleration (`RAPIDS-singlecell`) for integration and clustering steps. If `false` the value will be used in the configuration file for new projects created with the `scprocess newproj -c` command.

##### arvados

* `arv_instance` (optional): the name of the default Arvados instance. If specified it will be used in the configuration file for new projects created with the `scprocess newproj -c` command.

##### ref_txomes

Prebuilt human and mouse reference transcriptomes from 10x Genomics can be downloaded with {{scsetup}} by adding `tenx` to the `scprocess_setup.yaml` file. Valid values for names are `human_2024`, `mouse_2024`, `human_2020`, `mouse_2020`.

Names and specifications for custom references should be listed in the `custom` section of the `scprocess_setup.yaml` file. For each `custom` genome users have to provide the following parameters:

* `name`: name to be used for the reference
* `fasta`: path to FASTA file
* `gtf`: path to GTF file
* `mito_str`: regular expression used to identify genes in the mitochondrial genome (example for mouse: `"^mt-"`)

Optional parameters for both `tenx` and `custom` references are:

* `decoys`: whether to include decoy sequences when building the mapping index; default `true`.

Optional parameter for `tenx` references is:

* `rrnas`: whether or not ribosomal RNAs should be included in the reference. If not specified it defaults to `true` for all `tenx` genomes.

??? note "Impact of custom parameters for `tenx` genomes on `scsetup` runtime"

    When configuring `tenx` genomes with their default values, {{scsetup}} will download prebuilt indices for `simpleaf`. However, if the default parameters are modified (e.g., setting `rrnas` or `decoys` to `false`), {{scsetup}} will build the indices from scratch during execution, which will increase the runtime.

??? info "More about decoys"
    {{sc}} utilizes `simpleaf`, a lightweight mapping approach that, by default, maps sequenced fragments exclusively to the transcriptome. However, this can lead to incorrect mapping of reads that arise from unannotated genomic loci to the transcriptome. To mitigate this issue, the `decoys` parameter for `ref_txomes` is set to `true`. This option allows `simpleaf` to identify genomic regions with sequences similar to those in transcribed regions (decoys), thereby reducing the likelihood of false mappings. For further details, refer to Srivastava et al., 2019[@Srivastava2020-jb].


##### probe_sets

Probe set indices for 10x Flex data are prepared with {{scsetup}} by adding a `probe_sets` section to the `scprocess_setup.yaml` file. Each entry under `tenx` specifies a named probe set to download and index. Valid values for `name` are `human_v1`, `human_v2`, `mouse_v1`, and `mouse_v2`.


##### xgboost { #setup-xgboost }

Custom `XGBoost` classifiers trained with {{sc}} (using the `train_xgboost` functionality) can be registered with {{scsetup}} so they are available for cell type annotation in any project. To register a user-trained classifier, add an `xgboost` section to the `scprocess_setup.yaml` file. Each entry specifies a classifier to install. There are two ways to reference a trained classifier:

* use the `config` parameter to point to the project configuration file that was used to train the classifier. {{sc}} will automatically locate the training outputs.
* use the `src_dir` parameter to point directly to the directory containing the trained classifier. Also specify `name`, the identifier to use when annotating other projects.

For a local classifier with display metadata, for example:

```yaml
xgboost:
  - name: my_classifier
    src_dir: /path/to/project/output/my_project_train_xgboost
    xgboost_tag: xgboost_my_classifier
```

For installation from `src_dir`, `xgboost_tag` is the full filename prefix of the trained model; `name` is the identifier used after installation. Alternatively, supply `config: /path/to/config.yaml` (and an optional `name`) to locate the model from its training configuration. 

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

**Description**: Create a new `workflowr` project directory for {{scjoin}} outputs.

**Parameters**:

* `name` (positional): name of the new `workflowr` project directory.
* `-w`/`--where` (optional): path to the parent directory where the new project will be created; defaults to the current working directory.
* `-c`/`--config` (optional): generates a template join configuration YAML file.
* `-p`/`--projects` (optional): one or more paths to existing {{sc}} project configuration files, used to pre-fill the generated config's `projects:` block. Requires `-c`; has no effect otherwise.


## {{scknee}} { #scprocess-plotknee }

**Description**: Create an interactive barcode-rank plot.

**Parameters**:

* `sample`: `sample_id` corresponding to the barcode-rank curve.
* `-k`/`--kneefile`: path to the knee plot data file generated by {{sc}}, e.g. `output/[short_tag]_mapping/af_[sample_id]/knee_plot_data_[sample_id]_[date_stamp].csv.gz`. Exactly one of `--kneefile` and `--configfile` should be specified.
* `-c`/`--configfile`: path to configuration file used for running {{sc}}. Exactly one of `--kneefile` and `--configfile` should be specified.


## {{scrun}} { #scprocess-run }

**Description**: Run {{sc}}.

**Parameters**:

* `configfile` (positional): path to a configuration YAML file.
* `-n`/`--dry-run` (optional): list the analysis steps without executing them. Helpful for checking input files and parameters.
* `-E`/`--extraargs` (optional): list of additional arguments to pass to `Snakemake`. Refer to [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for a detailed explanation of available command-line options.
* `--unlock` (optional): unlock the directory if a previous run was interrupted.
* `-r`/`--rule` (optional): specifies which rule {{sc}} should run. Options are:
    + `all`: default; runs [Core pipeline steps](introduction.md#core-pipeline-steps); `cell_cycle` and `label_celltypes` have to be configured.
    + `mapping`: read alignment and quantification.
    + `ambient`: ambient RNA removal (optional) and cell calling.
    + `demux`: sample demultiplexing.
    + `qc`: cell-level QC filtering and doublet identification.
    + `hvg`: calculation of highly variable genes.
    + `cell_cycle`: cell cycle anotation with `trycicle`.
    + `integration`: dimensionality reduction with PCA, optional batch correction with `Harmony`, UMAP and clustering.
    + `marker_genes`: marker gene identification and optional gene set enrichment analysis.
    + `label_celltypes`: cell type annotation.
    + `train_xgboost`: `XGBoost` classifier training.
    + `zoom`: subclustering.
    + `shiny`: interactive Shiny app generation.
* `--zoom` (optional): with `-r zoom`, `-r shiny` or `-r train_xgboost`, select a zoom (subclustering) by name (e.g. `--zoom immune_cells`). With `-r zoom` or `-r shiny`, `all` selects every zoom defined in the config; `-r train_xgboost` trains one classifier at a time, so `--zoom all` is not accepted.


### configuration file

This is an example config file for {{scrun}} with all parameters and their default values/placeholders. Required parameters are highlighted:

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
      tenx_assay_type: poly_a
      ref_txome: # should only be defined if tenx_assay_type is poly_a
      probe_set: # should only be defined if tnex_assay_type is flex
      tenx_chemistry: none
      metadata_vars:
      show_arv_uuids: true
      custom_sample_params:
      exclude:
        sample_id:
        pool_id:
    multiplexing:
      demux_type: none
      fastq_dir:
      arv_uuids:
      feature_ref:
      hto_id_regex:
      demux_output:
      seurat_quantile: 0.99
    ambient:
      ambient_method: decontx_background
      cell_calling: barcodeRanks
      cb_version: v0.3.2
      cb_max_prop_kept: 0.9
      cb_learning_rate: 0.0001
      cb_posterior_batch_size: 128
      cb_empty_training_fraction: 0.2
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
      dbl_min_counts: 0
      dbl_dbr_sd:
      dbl_clusters:
      exclude_mito: true
    pb_empties:
      ambient_genes_logfc_thr: 0
      ambient_genes_fdr_thr: 0.01
    hvg:
      hvg_method: sample
      hvg_n_hvgs: 2000
      hvg_chunk_size: 2000
      hvg_exclude_ambient_genes: true
      hvg_exclude_from_file:
      hvg_metadata_split_var:
    integration:
      int_pca_method: bpcells
      int_use_gpu: true
      int_embedding: harmony
      int_cl_method: leiden
      int_theta: 0.1
      int_batch_var: sample_id
      int_n_dims: 50
      int_umap_n_iter: 2000
      int_dbl_res: 4
      int_dbl_cl_prop:
      int_dbl_threshold_method: adaptive
      int_sce_outs: false
      int_use_paga: true
      int_paga_cl_res: 2
      int_cell_cycle_harmonics: 2
      int_cell_cycle_ridge_lambda: 0.1
      int_cell_cycle_correction:
      int_res_ls: [0.1, 0.2, 0.5, 1, 2]
    cell_cycle:
      cyc_density_percentile:
      cyc_density_percentile_grid: [0, 0.1, 0.3, 1, 3, 10, 100]
      cyc_bandwidth_multiplier: 1.0
      cyc_kde_grid_size: 500
      cyc_seed: 20230308
      cyc_origin_manual:
    marker_genes:
      mkr_sel_res: 0.2
      mkr_min_cl_size: 100
      mkr_min_cells: 10
      mkr_dispersion_use_clusters: true
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
        label_specs_f:
        coarse_label_specs_f:
        save_cluster_names_file: false
    train_xgboost:
      annots_f:
      label_specs_f:
      coarse_label_specs_f:
      min_cells_per_label: 20
      min_cells_expressed: 10
      gene_exclude_re: "(lincRNA|lncRNA|pseudogene|antisense)"
      seed: 42
      use_gpu: false
      pass1_subsample: 0.632
      pass1_colsample_bytree: 0.1
      pass1_learning_rate: 0.1
      pass1_nrounds: 1000
      pass1_early_stopping: 10
      pass2_colsample_bytree: 0.5
      pass2_learning_rate: 0.05
      pass2_nrounds: 2000
      pass2_early_stopping: 10
      downsample_copies: 0
      downsample_min_umis: 300
      downsample_max_umis: 3000
      gain_threshold: 0.9
      min_genes: 100
      max_genes: 5000
      xgboost_tag:
      excluded_labels: []
      cluster_majority_vote_resolution:
      cluster_majority_vote_min_purity: 0.5
      test_group_var:
      test_groups:
      validation_group_var: batch
      validation_fraction: 0.2
      validation_balancing_groups: []
      sampling_balancing_groups: []
      train_cells_per_label: 2000
      validation_cells_per_label: 500
    zoom:
    resources:
      retries: 0
      n_run_mapping: 8
      n_run_mapping_flex: 8
      n_run_marker_genes: 8
    shiny:
      app_tag:
      app_title:
      email:
      keyword: cells
      default_gene:
      n_keep: 30000
      metadata_vars:
      metadata_labels:
      metadata_combns:
      home_md:
      annotation_csv:
      cluster_palette:
      metadata_palettes:
    ```

=== "placeholders"

    ```yaml hl_lines="1 2 3 4 5 6 7 8 9 10 11 12 13 14"
    project:
      proj_dir: /path/to/proj/directory
      fastq_dir: /path/to/directory/with/fastq/files
      arv_uuids: ["instance_name1", "instance_name2", "instance_name3"]
      arv_instance: instance_name
      full_tag: test_project
      short_tag: test
      your_name: Test McUser
      affiliation: where you work
      date_stamp: "2050-01-01"
      sample_metadata: /path/to/metadata.csv
      tenx_assay_type: poly_a
      ref_txome: human_2024
      probe_set: human_v1
      tenx_chemistry: 3v3
      metadata_vars: [var1, var2]
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
      arv_uuids: ["instance_name1", "instance_name2", "instance_name3"]
      feature_ref: /path/to/feature_ref.csv
      hto_id_regex: "^HTO"
      demux_output: /path/to/demux_output.csv
      seurat_quantile: 0.99
    ambient:
      ambient_method: cellbender
      cell_calling: barcodeRanks
      cb_version: v0.3.2
      cb_max_prop_kept: 0.9
      cb_learning_rate: 0.001
      cb_posterior_batch_size: 128
      cb_empty_training_fraction: 0.2
      cb_expected_cells: 10000
      cb_total_droplets_included: 20000
      cb_low_count_threshold: 5
    qc:
      qc_min_counts: 500
      qc_max_counts: 20000
      qc_min_feats: 300
      qc_max_feats: 6000
      qc_min_mito: 0
      qc_max_mito: 0.1
      qc_min_splice: 0
      qc_max_splice: 1
      qc_min_cells: 100
      dbl_min_feats: 100
      dbl_min_counts: 0
      dbl_dbr_sd: 1
      dbl_clusters: true
      exclude_mito: true
    pb_empties:
      ambient_genes_logfc_thr: 0
      ambient_genes_fdr_thr: 0.01
    hvg:
      hvg_method: sample
      hvg_n_hvgs: 2000
      hvg_chunk_size: 2000
      hvg_exclude_ambient_genes: true
      hvg_exclude_from_file: /path/to/file/with/genes/to/exclude
      hvg_metadata_split_var: var1
    integration:
      int_pca_method: bpcells
      int_use_gpu: true
      int_embedding: harmony
      int_cl_method: leiden
      int_theta: 0.1
      int_batch_var: sample_id
      int_n_dims: 50
      int_umap_n_iter: 2000
      int_dbl_res: 4
      int_dbl_cl_prop: 0.5
      int_dbl_threshold_method: adaptive
      int_sce_outs: false
      int_use_paga: true
      int_paga_cl_res: 2
      int_cell_cycle_harmonics: 2
      int_cell_cycle_ridge_lambda: 0.1
      int_cell_cycle_correction: within_sample
      int_res_ls: [0.1, 0.2, 0.5, 1, 2]
    cell_cycle:
      cyc_density_percentile: 1.3
      cyc_density_percentile_grid: [0, 0.1, 0.3, 1, 3, 10, 100]
      cyc_bandwidth_multiplier: 1.0
      cyc_kde_grid_size: 500
      cyc_seed: 20230308
      cyc_origin_manual: [0.1, 0.2]
    marker_genes:
      mkr_sel_res: 0.2
      mkr_min_cl_size: 100
      mkr_min_cells: 10
      mkr_dispersion_use_clusters: true
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
        label_specs_f: /path/to/label_specs.csv
        coarse_label_specs_f: /path/to/coarse_label_specs.csv
        save_cluster_names_file: true
    train_xgboost:
      annots_f: /path/to/annotations.csv.gz
      label_specs_f: /path/to/label_specs.csv
      coarse_label_specs_f: /path/to/coarse_label_specs.csv
      min_cells_per_label: 20
      min_cells_expressed: 10
      gene_exclude_re: "(lincRNA|lncRNA|pseudogene|antisense)"
      seed: 42
      use_gpu: false
      pass1_subsample: 0.632
      pass1_colsample_bytree: 0.1
      pass1_learning_rate: 0.1
      pass1_nrounds: 1000
      pass1_early_stopping: 10
      pass2_colsample_bytree: 0.5
      pass2_learning_rate: 0.05
      pass2_nrounds: 2000
      pass2_early_stopping: 10
      downsample_copies: 0
      downsample_min_umis: 300
      downsample_max_umis: 3000
      gain_threshold: 0.9
      min_genes: 100
      max_genes: 5000
      xgboost_tag: my_classifier
      excluded_labels: [unassigned]
      cluster_majority_vote_resolution: "RNA_snn_res.2"
      cluster_majority_vote_min_purity: 0.5
      test_group_var: sample_id
      test_groups: [sample1, sample2]
      validation_group_var: batch [**how can this be**?]
      validation_fraction: 0.2
      validation_balancing_groups: [treatment]
      sampling_balancing_groups: [treatment]
      train_cells_per_label: 2000
      validation_cells_per_label: 500
    zoom:
      - /path/to/cell_subset_1_zoom_params.yaml
      - /path/to/cell_subset_2_zoom_params.yaml
      - /path/to/cell_subset_3_zoom_params.yaml
    resources:
      retries: 3
      n_run_mapping: 8
      n_run_mapping_flex: 8
      n_run_marker_genes: 8
    shiny:
      app_tag: test_project
      app_title: My Project
      email: me@example.com
      keyword: cells
      default_gene: Mog
      n_keep: 20000
      metadata_vars: [var1, var2]
      metadata_labels:
        brainregion: "brain region"
        condition: "treatment group"
      metadata_combns: [[var1, var2]]
      home_md: /path/to/home.md
      annotation_csv: /path/to/annotation.csv
      cluster_palette: Set1
      metadata_palettes:
        condition: Dark2
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
* `tenx_assay_type`: assay type. Options are `poly_a` and `flex`. When set to `flex`, `probe_set` is required and `ref_txome` must not be specified.
* `ref_txome`: must match one of the values in the `ref_txome` column of `index_parameters.csv` (created by {{scsetup}}). Required for polyA data; must not be specified for Flex data (use `probe_set` instead).
* `probe_set`: probe set to use for 10x Flex data. Required for Flex data. Valid values are `human_v1`, `human_v2`, `mouse_v1`, and `mouse_v2`. Must not be specified for polyA data (use `ref_txome` instead).

#### Optional parameters

##### project

* `tenx_chemistry`: 10x assay configuration. Accepted values are `3v2`, `3v3`, `3v4`, `5v1`, `5v2`, `5v3`, `multiome`, `flexv1` and `flexv2`. `multiome` refers only to gene expression data generated with the 10x multiome kit (ATACseq data is not supported). For Flex data, chemistry specification is redundant as chemistry is derived automatically based on `probe_set`. When not specified for polyA data, chemistry is auto-detected.
* `metadata_vars`: a list of column names in the `sample_metadata` file to be used for visualizing the distribution of cell annotations across identified clusters and regions of the low-dimensional embedding.
* `show_arv_uuids`: whether to display Arvados UUIDs (`arv_uuids`) in the configuration file details box on the index page. If `false`, UUIDs are replaced with "not shown".
* `exclude`: list of all samples that should be excluded from the analysis. Samples can be listed under `pool_id` (if multiplexed) or `sample_id`.
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
      cb_posterior_batch_size: 128 # only applicable if cb_version is v0.3.2
sample_id:
  sample_1:
    qc:
      qc_min_counts: 100
      dbl_dbr_sd: 1
```

##### multiplexing

* `demux_type`: multiplexing approach used. Options are:
    + `none` if experiment is not multiplexed;
    + `hto` if hto-based demultiplexing of samples should be performed with {{sc}};
    + `ocm` if on-chip multiplexing was used. Requires an `ocm_id` column in sample metadata with values `OB1`-`OB4`.
    + `custom` if demultiplexing results will be used as input to {{sc}}; or
    + `flex` for 10x Flex data where samples are pooled within a library and demultiplexed using probe barcodes. Requires `tenx_assay_type: flex` in the `project` section.
* `fastq_dir`: path to directory containing HTO FASTQ files. Should be absolute or relative to `proj_dir`. If `demux_type` is `hto`, exactly one of `fastq_dir` and `arv_uuids` should be specified.
* `arv_uuids`: list of Arvados UUIDs where HTO fastq files are located. Expects `arv_instance` to be defined. If `demux_type` is `hto`, exactly one of `fastq_dir` and `arv_uuids` should be specified.
* `feature_ref`: path to CSV file with columns `hto_id` and `sequence`. Required if `demux_type` is `hto`.
* `hto_id_regex`: regular expression for the `hto_id` values of the HTOs in `feature_ref`. Rows that do not match are condiered to correspond to other antibody-derived tags used in the experiment. Used if `demux_type` is `hto`.
* `demux_output`: path to CSV file with columns `pool_id`, `sample_id`, `cell_id`. Optional column `class` can be added with values `doublet`, `singlet` or `negative`. Required if `demux_type` is `custom`.
* `seurat_quantile`: equivalent to the `positive.quantile` argument of the `Seurat::HTODemux` function (see [Seurat documentation](https://satijalab.org/seurat/reference/htodemux) for more details). Used if `demux_type` is `hto`.

##### ambient

* `ambient_method`: method for ambient RNA removal. Options are:
    + `decontx_background`: runs [`DecontX`](https://bioconductor.org/packages/release/bioc/html/celda.html) with an estimated ambient RNA profile as input.
    + `decontx_cluster`: runs `DecontX` with a cluster-based background model.
    + `cellbender`: uses [CellBender](https://cellbender.readthedocs.io) for ambient RNA removal.
    + `none`: skips ambient RNA removal;
* `cell_calling`: method for cell calling when `ambient_method` is `none`, `decontx_background`, or `decontx_cluster`. Options are `barcodeRanks` and `emptyDrops`.
* `cb_version`: CellBender version to use with `ambient_method: cellbender`; Options are `v0.3.2`, `v0.3.0` and `v0.2.0`.
* `cb_max_prop_kept`: maximum proportion of droplets, relative to `--total-droplets-included`, that `cellbender` can call as cells — samples are excluded if `cellbender` calls more than this proportion of `--total-droplets-included` droplets as cells. Applicable only if `ambient_method` is `cellbender`. For more information about the `--total-droplets-included` parameter see [Cellbender documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_learning_rate`: sets the `--learning-rate` `CellBender` parameter to the specified value; applicable only if `ambient_method` is `cellbender`. For more information about this parameter see [Cellbender documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_empty_training_fraction`: sets the `--empty-drop-training-fraction` `CellBender` parameter to the specified value; applicable only if `ambient_method` is `cellbender`. Setting this to a lower value (e.g. 0.1 or 0.05) can help if `CellBender` jobs are failing on samples with very few cells. For more information about this parameter see [Cellbender documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_posterior_batch_size`: value of the `--posterior-batch-size` parameter; applicable only if `ambient_method` is `cellbender` and `cb_version` is `v0.3.2`. For more information about this parameter see [Cellbender documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_expected_cells`: forces the `--expected-cells` `Cellbender` parameter to be consistent across all samples; applicable only if `ambient_method` is `cellbender`. For more information about this parameter see [Cellbender documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_total_droplets_included`: forces the `--total-droplets-included` `Cellbender` parameter to be consistent across all samples; applicable only if `ambient_method` is `cellbender`. For more information about this parameter see [Cellbender documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_low_count_threshold`: forces the `--low-count-threshold` `CellBender` parameter to be consistent across all samples; applicable only if `ambient_method` is `cellbender`. For more information about this parameter see [Cellbender documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).


##### qc

* `qc_min_counts`: minimum number of UMIs per cell required to retain a cell.
* `qc_max_counts`: maximum number of UMIs per cell required to retain a cell. No upper limit if unset.
* `qc_min_feats`: minimum number of detected features per cell required to retain a cell.
* `qc_max_feats`: maximum detected features per cell required to retain a cell. No upper limit if unset.
* `qc_min_mito`: minimum proportion of mitochondrial reads required to retain a cell.
* `qc_max_mito`: maximum proportion of mitochondrial reads required to retain a cell.
* `qc_min_splice`: minimum proportion of spliced reads required to retain a cell.
* `qc_max_splice`: maximum proportion of spliced reads required to retain a cell.
* `qc_min_cells`: minimum cells remaining in a sample after QC to retain that sample.
* `exclude_mito`: boolean; whether to exclude mitochondrial genes or not.
* `dbl_clusters`: whether artificial doublets are generated using cluster labels. If unset, uses the `scDblFinder` default. `false` disables cluster-based generation; `true` first clusters the cells and generates artificial doublets between clusters.
* `dbl_dbr_sd`: how strongly the expected doublet rate constrains the calling threshold. If unset, uses the `scDblFinder` default. Setting `1` removes the expected-rate penalty, leaving the threshold determined by classification of artificial doublets.
* `dbl_min_feats`: minimum number of features required for each barcode to be included in `scDblFinder` calculations.
* `dbl_min_counts`: minimum number of UMIs required for each barcode to be included in `scDblFinder` calculations; `0` disables this filter.


##### pb_empties

* `ambient_genes_logfc_thr`: log-fold change (logFC) threshold used to filter the results of the `edgeR` differential expression test comparing empty droplets to cells.
* `ambient_genes_fdr_thr`: false discovery rate (FDR) threshold used to filter the results of the `edgeR` differential expression test comparing empty droplets to cells.

##### hvg

* `hvg_method`: options:
    + `sample`: calculate highly variable genes per sample, then calculate combined ranking across samples;
    + `all`: calculate highly variable genes across all cells in the dataset; and
    + `groups`: calculate highly variable genes for each sample group then calculate combined ranking across groups.
* `hvg_metadata_split_var`: if `hvg_method` is `groups`, which variable in `sample_metadata` should be used to define sample groups.
* `hvg_n_hvgs`: number of HVGs to use for PCA.
* `hvg_exclude_ambient_genes`: if `true`, genes enriched in empty droplets relative to cells will be excluded from highly variable genes selection.
* `hvg_exclude_from_file`: path to CSV file with genes to be excluded from HVGs. Should be absolute or relative to `proj_dir`. File should contain one column, named either `gene_id` or `symbol`. Values in the column should all be present in reference genome.
* `hvg_chunk_size`: number of genes to use for each chunked matrix.

##### cell_cycle

* `cyc_density_percentile` (required): percentile (0–100) of cell density used as the cutoff for full weight. Cells at or below this density receive full weight, while cells in denser regions are down-weighted. At `0`, weights are proportional to inverse density; at `100`, every cell receives equal weight.
* `cyc_origin_manual`: manually selected cell cycle origin.
* `cyc_density_percentile_grid`: alternative density percentiles shown in the report.
* `cyc_bandwidth_multiplier`: positive multiplier controlling how broadly neighbouring cells are grouped when estimating the density of cells the projection.
* `cyc_kde_grid_size`: positive integer specifying the number of grid bins per axis used to evaluate the 2D kernel density estimate across the projection (an N × N grid).
* `cyc_seed`: random seed for generating diagnostic plots.[**Change this default seed**?]

##### integration

* `int_embedding`: use `harmony` to correct for the chosen batch variable, or `pca` to skip batch correction.
* `int_batch_var`: the variable identifying batches for `Harmony`. If `demux_type` is set to either `hto` or `custom`, then `pool_id` is an alternative option.
* `int_theta`: theta parameter for `Harmony` integration, controlling batch mixing.
* `int_n_dims`: number of principal components to use for downstream analysis.
* `int_res_ls`: list of resolution values to use for clustering.
* `int_dbl_res`: clustering resolution for identification of additional doublets.
* `int_dbl_threshold_method`: method for setting the cluster-level doublet threshold - clusters exceeding this proportion of doublets are removed. `adaptive` estimates the cutoff automatically; `manual` uses a user-defined cutoff (`int_dbl_cl_prop`).
* `int_dbl_cl_prop`: cutoff for the `manual` method. Required when `int_dbl_threshold_method` is `manual`. Should be omitted if `int_dbl_threshold_method` is `adaptive`.
* `int_umap_n_iter`: number of steps used to optimize the UMAP layout.
* `int_use_paga`: if `true`, enable Partition-based graph abstraction (PAGA) and initialize a UMAP using PAGA results.
* `int_paga_cl_res`: clustering resolution for PAGA analysis. Must be a value listed in `int_res_ls`. Only used when `int_use_paga` is `true`.
* `int_pca_method`: `bpcells` streams normalized, scaled expression data from disk for PCA. `in_memory` loads expression matrix into memory and computes PCA via `Scanpy` (CPU) or `RAPIDS-singlecell` (GPU, if `int_use_gpu` is set).
* `int_use_gpu`: boolean; whether to use GPU acceleration (`RAPIDS-singlecell`). If GPU is not available, `Scanpy` will be used.
* `int_sce_outs`: if `true` H5AD outputs will be converted to `SingleCellExperiment` objects and stored as RDS files.
* `int_cell_cycle_correction`: (optional; omitted = no correction) removes cell-cycle-phase variation from the PCA input using each cell's `tricycle_theta` score. Both options remove phase variation *within* each sample. `within_sample` leaves each sample's average expression unchanged, so cycling-composition differences between samples remain for `Harmony` to resolve. `global` also removes the part of that between-sample average difference attributable to differing cell-cycle composition — use it only if that dominates the leading PCs.
* `int_cell_cycle_harmonics`: number of sine/cosine pairs used to fit expression against cell-cycle phase; `1` fits a single wave around the cycle, `2` adds a second harmonic so the fit can have two peaks and troughs per cycle, closer to how many cell-cycle genes actually behave.
* `int_cell_cycle_ridge_lambda`: ridge penalty shrinking the fitted cell-cycle correction's coefficients (see `int_cell_cycle_harmonics`) toward zero; larger values give a weaker, more conservative correction, and also protect against a degenerate fit when a sample/group has too few cells to estimate it from reliably. `0` disables shrinkage.


##### marker_genes

* `mkr_sel_res`: clustering resolution used for marker analysis; Must match one of the elements in `int_res_ls`.
* `mkr_min_cl_size`: minimum number of cells required in a cluster to calculate marker genes for that cluster.
* `mkr_min_cells`: minimum number of cells required in a pseudobulk sample to include it in marker gene calculations.
* `mkr_min_cpm_mkr`: minimum expression in counts per million (CPM) for selecting genes shown in marker plots.
* `mkr_not_ok_re`: pattern matching gene types to omit from marker-expression plots.
* `mkr_do_gsea`: boolean specifying whether Gene Set Enrichment Analysis (GSEA) should be performed on marker genes.
* `mkr_gsea_cut`: false discovery rate (FDR) cutoff for GSEA.
* `mkr_min_cpm_go`: minimum CPM in a cluster required for a gene to be used in GSEA.
* `mkr_max_zero_p`: maximum proportion of pseudobulk samples for a cluster that can have zero counts for a gene to be used in GSEA.
* `mkr_gsea_var`: the statistic used for ranking genes in GSEA. Choices are `z_score` (z-score based on signed `log10(FDR)`) or `logFC` (log fold change).
* `mkr_custom_genesets`: a list of custom marker gene sets, each defined by a unique name and associated file path.
    + `name`: a string representing the name of the marker gene set
    + `file`: path to CSV file containing a list of genes in the marker gene set. Must contain column `label` (marker gene category), and `symbol` and/or `ensembl_id`. If not speficied `scprocess` will look for file `$SCPROCESS_DATA_DIR/marker_genes/{name}.csv`
* `mkr_dispersion_use_clusters`: boolean specifying whether to include cluster membership when estimating variation between samples. [**Where is this parameter relevant??**]


##### label_celltypes

* `labeller`: specifies the method to annotate cell types; options include:
    + `celltypist` uses one of the models available in [`CellTypist`](https://www.celltypist.org/models) for annotation.
    + `scprocess`: uses an `XGBoost` classifier for cell type annotation.
* `model`: determines the model to be used based on the selected `labeller`. For list of all available `CellTypist` models see `$SCPROCESS_DATA_DIR/celltypist/celltypist_models.csv`. For a list of all available `scprocess` classifiers see `$SCPROCESS_DATA_DIR/xgboost/xgboost_models.csv`.
* `hi_res_cl`: name of a column containing high-resolution clustering results. It must follow the pattern `"RNA_snn_res.n"` where `n` should be replaced with one of the values in `int_res_ls`.
* `min_cl_prop`: minimum fraction of a cluster sharing its most frequent predicted label for that label to be assigned to the whole cluster.
* `min_cl_size`: minimum number of cells in a cluster required for that cluster to be labeled.
* `label_specs_f`: CSV file with one row per label, listed in the desired display order. Required column: `label`. Optional columns: `colour` and `coarse_label`. Should be absolute or relative to `proj_dir`. If omitted, {{sc}} looks for `$SCPROCESS_DATA_DIR/celltypist/{model}_label_specs.csv` (if `labeller` is `celltypist`) or `$SCPROCESS_DATA_DIR/xgboost/{model}/{model}_label_specs.csv` (if `labeller` is `scprocess`);
* `coarse_label_specs_f`: CSV file with one row per coarse group, listed in the desired display order. Required column: `coarse_label`. Optional column: `coarse_colour`. Requires a coarse mapping in `label_specs_f`. Resolved the same way as `label_specs_f`: an explicit value takes priority, otherwise {{sc}} looks for the matching `{model}_coarse_label_specs.csv` under `$SCPROCESS_DATA_DIR`.
* `save_cluster_names_file`: boolean specifying whether to generate a `cluster_names_for_shiny_*.csv` file mapping clusters (at `mkr_sel_res` resolution) to predicted cell type names. This file can be used as the `annotation_csv` in the shiny config.


##### train_xgboost

??? info "Splitting the input dataset into training, validation, and test sets"

    Prior to training the classifier, cells are partitioned into training, validation, and test sets at the group level. All cells belongng to a given group are assigned to the same set.

    By default, groups correspond to values of `int_batch_var` (`sample_id` or `pool_id`). Groups can also be defined by any other metadata variable.

    - **Test set:** If `test_group_var` is specified, all samples corresponding to `test_groups` form the test set. These cells are fully held out from training and are reserved for diagnostic plots and evaluation.
    - **Validation set:** If `validation_group_var` is specified, a fraction of its unique values (`validation_fraction`) is randomly selected, and all samples with those values form the validation set. Unlike test cells, validation cells influence training, for example by determining early stopping during boosting.
    - **Training set:** All remaining cells are used for training, with up to `train_cells_per_label` cells sampled per class.

* `annots_f` (required): CSV file with a unique `cell_id` column and a `label` column containing the training labels. It may also include other metadata columns, but any values shared with the `sample_metadata` must match. Should be absolute or relative to `proj_dir`. [**sould we allow other metadata variables in this file that are not in sample_metadata? - it just makes things more complicated**]
* `label_specs_f` (optional): CSV file with one row per label, listed in the desired display order. Required column: `label`. Optional columns: `colour` and `coarse_label`.
* `coarse_label_specs_f`: CSV file with one row per coarse group, listed in the desired display order. Required column: `coarse_label`. Optional column: `coarse_colour`. Requires a coarse mapping in `label_specs_f`.
* `cluster_majority_vote_resolution`: a clustering resolution at which each cluster's most frequent training label is assigned to all cells in that cluster. Omit to keep the original labels.
* `cluster_majority_vote_min_purity`: minimum fraction of labelled training cells that must share the most frequent label before it is assigned to the cluster. [**Only training cells are used for the vote, so validation and test labels do not affect the decision.**]
* `excluded_labels`: labels to exclude from training.
* `min_cells_per_label`: minimum number of cells required to keep a label for training. Labels with fewer cells are excluded.
* `test_group_var`: metadata column that defines the values to hold out completely from training and validation. Must be supplied together with `test_groups`.
* `test_groups`: values of `test_group_var`.
* `validation_group_var`: metadata column that defines the groups kept together when splitting data for validation. [**How can the default be `batch`?**]
* `validation_fraction`: fraction of groups assigned to the validation set.
* `validation_balancing_groups`: metadata columns (other than `validation_group_var`) representing sample groups that should appear in both the validation and training sets. If multiple columns are given, every combination of their values is identified, and each combination is ensured to contribute to both sets.
* `train_cells_per_label`: maximum number of cells selected per label after groups are assigned.
* `validation_cells_per_label`: maximum number of cells selected per label after groups are assigned.
* `sampling_balancing_groups`: if a label has more cells than `train_cells_per_label`, cell selection is balanced within that label across these metadata variables.
* `min_cells_expressed`: minimum number of training cells in which a gene must be detected.
* `gene_exclude_re`: regular expression matching gene types to exclude.
* `gain_threshold`: fraction of total gain from the first fit that the selected genes must account for.
* `min_genes`: lower limit on the number of selected genes.
* `max_genes`: upper limit on the number of selected genes.
* `seed`: random seed for reproducibility.
* `use_gpu`: whether to train on a GPU.
* `xgboost_tag`: name used in model filenames. If omitted, the file prefix is `xgboost_{full_tag}` for a project or `xgboost_{join.name}` for a join.
* `pass1_subsample`: fraction of training cells sampled for each boosting round in the first fit.
* `pass1_colsample_bytree`: fraction of genes considered for each tree in the first fit.
* `pass1_learning_rate`: contribution of each added tree to the predictions in the first fit.
* `pass1_nrounds`: maximum number of boosting rounds in the first fit.
* `pass1_early_stopping`: stop the first fit after this many rounds without improvement in validation performance.
* `pass2_colsample_bytree`: fraction of genes considered for each tree in the final fit. Only the genes selected after the first fit are considered.
* `pass2_learning_rate`: contribution of each added tree to the predictions in the final fit.
* `pass2_nrounds`: maximum number of boosting rounds in the final fit.
* `pass2_early_stopping`: stop the final fit after this many rounds without improvement in validation performance.
[**these need further review**]
* `downsample_copies`: shallow copies per selected training and validation cell; `0` disables this. `1` is a good starting point.
* `downsample_min_umis`: start of range of copy depths in UMIs.
* `downsample_max_umis`: end of range of copy depths in UMIs.

##### shiny

* `app_tag`: app identifier used for directory and file names. 
* `app_title`: title displayed in the Shiny app header.
* `email`: contact email shown in the app footer.
* `keyword`: short word used in plot axis labels and descriptions (e.g. `"cells"`, `"nuclei"`).
* `default_gene`: gene symbol displayed by default in the "Explore Genes" tab.
* `n_keep`: number of cells retained in the subsampled UMAP shown in the app.
* `metadata_vars`: metadata variables to include in the app's plots. If not specified, project-level `metadata_vars` are used.
* `metadata_labels`: display labels for metadata variables. Keys are column names, values are display labels. Variables not listed keep their column name as the label.
* `metadata_combns`: list of metadata variable pairs to display as combined groupings. Each element should be a two-element list of variable names from `metadata_vars`.
* `home_md`: path to a Markdown file used as the landing page content. Should be absolute or relative to `proj_dir`.[**Is that true for join**?]
* `annotation_csv`: path to a CSV file with columns `cluster`, `cluster_name`, and optionally `colour`, defining display names, order, and colours for clusters. Should be absolute or relative to `proj_dir`.
* `cluster_palette`: name of colour palette applied to clusters when `annotation_csv` is not provided. Accepts any name from the Supported colour palettes list.
* `metadata_palettes`: per-variable colour palette configuration. Each key is a metadata variable name; the value can be:
    - a palette name (string)
    - an explicit colour list (array)
    - an object with optional `palette`, `colours`, and `values` keys (the `values` list controls level ordering; defaults to frequency order)[**This needs examples**]

??? info "Supported colour palettes"

    The following palette names are accepted anywhere a `cluster_palette` or `metadata_palettes` entry expects a name:

    | Source | Names |
    |--------|-------|
    | **scprocess built-in** | `nice_cols` (42 colours) |
    | **MetBrewer** | `Archambault`, `Austria`, `Benedictus`, `Cassatt1`, `Cassatt2`, `Cross`, `Degas`, `Demuth`, `Derain`, `Egypt`, `Gauguin`, `Greek`, `Hiroshige`, `Hokusai1`, `Hokusai2`, `Hokusai3`, `Homer1`, `Homer2`, `Ingres`, `Isfahan1`, `Isfahan2`, `Java`, `Johnson`, `Juarez`, `Kandinsky`, `Klimt`, `Lakota`, `Manet`, `Monet`, `Moreau`, `Morgenstern`, `Nattier`, `Navajo`, `NewKingdom`, `Nizami`, `OKeeffe1`, `OKeeffe2`, `Paquin`, `Peru1`, `Peru2`, `Pillement`, `Pissaro`, `Redon`, `Renoir`, `Signac`, `Tam`, `Tara`, `Thomas`, `Tiepolo`, `Troy`, `Tsimshian`, `VanGogh1`, `VanGogh2`, `VanGogh3`, `Veronese`, `Wissing` |
    | **RColorBrewer** | Qualitative: `Accent`, `Dark2`, `Paired`, `Pastel1`, `Pastel2`, `Set1`, `Set2`, `Set3`. Sequential/diverging: `Blues`, `BuGn`, `BuPu`, `GnBu`, `Greens`, `Greys`, `Oranges`, `OrRd`, `PuBu`, `PuBuGn`, `PuRd`, `Purples`, `RdPu`, `Reds`, `YlGn`, `YlGnBu`, `YlOrBr`, `YlOrRd`, `BrBG`, `PiYG`, `PRGn`, `PuOr`, `RdBu`, `RdGy`, `RdYlBu`, `RdYlGn`, `Spectral`. |
    | **ggsci** | `npg`, `aaas`, `nejm`, `lancet`, `jama`, `jco`, `ucscgb`, `d3` / `d3_20`, `d3_10`, `d3_20b`, `d3_20c`, `igv`, `cosmic`, `simpsons`, `rickandmorty`, `futurama`, `tron`, `startrek`, `uchicago`, `frontiers`, `flatui`, `bootstrap` |


##### zoom

In this section, users can provide multiple YAML files, each specifying parameters for repeating certain steps of {{sc}} on a subset of cells. The zoom configs can optionally contain `qc`, `hvg`, `integration`, `marker_genes`, `train_xgboost`, and `shiny` sections with same parameters as main project-level config file, with the following exceptions:

* `qc`: `dbl_min_feats`, `dbl_min_counts`, `dbl_dbr_sd`, `dbl_clusters`, `exclude_mito`
* `integration`: `int_batch_var`, `int_dbl_res`, `int_dbl_cl_prop`, `int_dbl_threshold_method`, `int_sce_outs`
* `shiny`: `app_tag`

Additional parameters include:

* `name`: name of cell subset to be analysed.
* `labels_source`: specifies how a cell subset is defined (required). Options include:
    - `scprocess`: labels assigned with the `XGBoost` classifier (using rule `label_celltypes`)
    - `celltypist`: labels assigned with `CellTypist`(using the rule `label_celltypes`)
    - `clusters`: labels based on clustering results obtained with {{sc}}
    - `custom`: user-defined cell type annotations
* `model`: required if `labels_source` is `scprocess` or `celltypist`.
* `sel_labels`: a list of all labels that define cell types/clusters to be included in subclustering (required).
* `labels_col`: name of the column containing the labels used to select cells (required).
* `save_subset_sces`: whether to create `RDS` files with `SingleCellExperiment` objects containing cells that have been assigned one of the values in `sel_labels`.
* `save_subset_anndata`: whether to create `H5AD` files containing cells that have been assigned one of the values in `sel_labels`.
* `custom_labels_f`: required if `labels_source` is set to `custom`; path to CSV file with columns `sample_id`, `cell_id` and `label`. Should be absolute or relative to `project_dir`. 
* `save_cluster_names_file`: boolean specifying whether to generate a `cluster_names_for_shiny_*.csv` file mapping this zoom's clusters (at `mkr_sel_res` resolution) to predicted cell type names, for use as the `annotation_csv` in the `shiny` config.
* `exclude`: samples or pools to drop from just the zoom's cell subset. 

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
  exclude:
    sample_id:
      - sample1
      - sample2
    pool_id:
      - pool1
      - pool2
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

This section allows users to adjust the resource requirements for specific `Snakemake` rules in {{scrun}}. This is especially useful when a step/rule fails on a cluster due to insufficient memory or runtime limits. By specifying the parameters below, users can fine-tune these settings for their pipeline:

* `gb_[rule_name]`: specifies the maximum memory (in GB) requested for running a specific rule. `rule_name` should be replaced with an {{sc}} rule name. This value applies for the entire job, not per thread.
* `mins_[rule_name]`: specifies the maximum runtime (in minutes) requested for running a specific rule. `rule_name` should be replace with an {{sc}} rule name.

Additional parameters include:

* `retries`: number of times to retry running a specific rule in {{scrun}} if it fails. For each attempt the initial memory requested for the rule is multiplied by `1.5**(attempt - 1)`. Useful for when {{scrun}} is ran on a [cluster](setup.md#cluster-setup).
* `n_run_mapping`: number of threads requested for running the mapping step if `tenx_assay_type` is `poly_a`.
* `n_run_mapping_flex`: number of threads requested for running the mapping step if `tenx_assay_type` is `flex`.
* `n_run_marker_genes`: number of threads requested for rule `run_marker_genes`.


??? note "Detailed information about resource parameters"

    * `gb_run_mapping_flex`: maximum memory required (in GB) for rule `run_mapping_flex`.
    * `gb_run_mapping`: maximum memory required (in GB) for rule `run_mapping`.
    * `gb_save_alevin_to_h5`: maximum memory required (in GB) for rule `save_alevin_to_h5`.
    * `gb_build_hto_index`: maximum memory required (in GB) for rule `build_hto_index`.
    * `gb_run_mapping_hto`: maximum memory required (in GB) for rule `run_mapping_hto`.
    * `gb_save_alevin_hto_to_h5`: maximum memory required (in GB) for rule `save_alevin_hto_to_h5`.
    * `gb_save_alevin_adt_to_h5`: maximum memory required (in GB) for rule `save_alevin_adt_to_h5`.
    * `gb_make_hto_sce_objects`: maximum memory required (in GB) for rule `make_hto_sce_objects`.
    * `gb_run_cellbender`: maximum memory required (in GB) for rule `run_cellbender`.
    * `gb_run_decontx`: maximum memory required (in GB) for rule `run_decontx`.
    * `gb_run_cell_calling`: maximum memory required (in GB) for rule `run_cell_calling`.
    * `gb_get_barcode_qc_metrics`: maximum memory required (in GB) for rule `get_barcode_qc_metrics`.
    * `gb_check_demux_ids`: maximum memory required (in GB) for rule `check_demux_ids`.
    * `gb_run_doublet_id_one_run`: maximum memory required (in GB) for rule `run_doublet_id_one_run`.
    * `gb_filter_qc_one_run`: maximum memory required (in GB) for rule `filter_qc_one_run`.
    * `gb_merge_doublet_diagnostics`: maximum memory required (in GB) for rule `merge_doublet_diagnostics`.
    * `gb_merge_qc`: maximum memory required (in GB) for rule `merge_qc`.
    * `gb_merge_rowdata`: maximum memory required (in GB) for rule `merge_rowdata`.
    * `gb_get_qc_sample_statistics`: maximum memory required (in GB) for rule `get_qc_sample_statistics`.
    * `gb_make_one_pb_empty`: maximum memory required (in GB) for rule `make_one_pb_empty`.
    * `gb_merge_pb_empty`: maximum memory required (in GB) for rule `merge_pb_empty`.
    * `gb_calculate_ambient_genes`: maximum memory required (in GB) for rule `calculate_ambient_genes`.
    * `gb_make_one_pb_cells`: maximum memory required (in GB) for rule `make_one_pb_cells`.
    * `gb_merge_pb_cells`: maximum memory required (in GB) for rule `merge_pb_cells`.
    * `gb_make_tmp_csr_matrix`: maximum memory required (in GB) for rule `make_tmp_csr_matrix`.
    * `gb_get_stats_for_std_variance_for_sample`: maximum memory required (in GB) for rule `get_stats_for_std_variance_for_sample`.
    * `gb_get_mean_var_for_group`: maximum memory required (in GB) for rule `get_mean_var_for_group`.
    * `gb_get_estimated_variances`: maximum memory required (in GB) for rule `get_estimated_variances`.
    * `gb_get_stats_for_std_variance_for_group`: maximum memory required (in GB) for rule `get_stats_for_std_variance_for_group`.
    * `gb_get_highly_variable_genes`: maximum memory required (in GB) for rule `get_highly_variable_genes`.
    * `gb_create_hvg_matrix`: maximum memory required (in GB) for rule `create_hvg_matrix`.
    * `gb_create_doublets_hvg_matrix`: maximum memory required (in GB) for rule `create_doublets_hvg_matrix`.
    * `gb_calculate_tricycle_projection_one_run`: maximum memory required (in GB) for rule `calculate_tricycle_projection_one_run`.
    * `gb_aggregate_tricycle`: maximum memory required (in GB) for rule `aggregate_tricycle`.
    * `gb_estimate_tricycle_reference_ellipse`: maximum memory required (in GB) for rule `estimate_tricycle_reference_ellipse`.
    * `gb_run_integration`: maximum memory required (in GB) for rule `run_integration`.
    * `gb_make_clean_h5ads`: maximum memory required (in GB) for rule `make_clean_h5ads`.
    * `gb_convert_h5ad_to_sce`: maximum memory required (in GB) for rule `convert_h5ad_to_sce`.
    * `gb_make_marker_pseudobulks`: maximum memory required (in GB) for rule `make_marker_pseudobulks`.
    * `gb_prepare_marker_pseudobulks`: maximum memory required (in GB) for rule `prepare_marker_pseudobulks`.
    * `gb_calc_marker_hvgs`: maximum memory required (in GB) for rule `calc_marker_hvgs`.
    * `gb_run_marker_genes`: maximum memory required (in GB) for rule `run_marker_genes`.
    * `gb_prepare_marker_plot_data`: maximum memory required (in GB) for rule `prepare_marker_plot_data`.
    * `gb_run_fgsea`: maximum memory required (in GB) for rule `run_fgsea`.
    * `gb_render_html_mapping`: maximum memory required (in GB) for rule `render_html_mapping`.
    * `gb_render_html_multiplexing`: maximum memory required (in GB) for rule `render_html_multiplexing`.
    * `gb_render_html_ambient`: maximum memory required (in GB) for rule `render_html_ambient`.
    * `gb_render_html_qc`: maximum memory required (in GB) for rule `render_html_qc`.
    * `gb_render_html_hvgs`: maximum memory required (in GB) for rule `render_html_hvgs`.
    * `gb_render_html_integration`: maximum memory required (in GB) for rule `render_html_integration`.
    * `gb_render_html_marker_genes`: maximum memory required (in GB) for rule `render_html_marker_genes`.
    * `gb_render_html_label_celltypes`: maximum memory required (in GB) for rule `render_html_label_celltypes`.
    * `gb_render_html_train_xgboost`: maximum memory required (in GB) for rule `render_html_train_xgboost`.
    * `gb_run_celltypist`: maximum memory required (in GB) for rule `run_celltypist`.
    * `gb_run_scprocess_labeller`: maximum memory required (in GB) for rule `run_scprocess_labeller`.
    * `gb_merge_labels`: maximum memory required (in GB) for rule `merge_labels`.
    * `gb_train_xgboost_prepare`: maximum memory required (in GB) for rule `train_xgboost_prepare`.
    * `gb_train_xgboost_train`: maximum memory required (in GB) for rule `train_xgboost_train`.
    * `gb_train_xgboost_predict`: maximum memory required (in GB) for rule `train_xgboost_predict`.
    * `gb_train_xgboost_combine`: maximum memory required (in GB) for rule `train_xgboost_combine`.
    * `gb_train_xgboost_depth_evaluation`: maximum memory required (in GB) for rule `train_xgboost_depth_evaluation`.
    * `gb_build_shiny_data`: maximum memory required (in GB) for rule `build_shiny_data`. The pre-split key `gb_build_shiny_app` is also accepted.
    * `gb_build_zoom_shiny_data`: maximum memory required (in GB) for rule `build_zoom_shiny_data`. The pre-split key `gb_build_zoom_shiny_app` is also accepted.
    * `gb_zoom_check_clusters`: maximum memory required (in GB) for rule `zoom_check_clusters`.
    * `gb_zoom_filter_cells_qc`: maximum memory required (in GB) for rule `zoom_filter_cells_qc`.
    * `gb_zoom_aggregate_labels`: maximum memory required (in GB) for rule `zoom_aggregate_labels`.
    * `gb_zoom_get_sample_statistics`: maximum memory required (in GB) for rule `zoom_get_sample_statistics`.
    * `gb_zoom_make_one_pb_cells`: maximum memory required (in GB) for rule `zoom_make_one_pb_cells`.
    * `gb_zoom_merge_pb_cells`: maximum memory required (in GB) for rule `zoom_merge_pb_cells`.
    * `gb_zoom_calculate_ambient_genes`: maximum memory required (in GB) for rule `zoom_calculate_ambient_genes`.
    * `gb_zoom_make_tmp_csr_matrix`: maximum memory required (in GB) for rule `zoom_make_tmp_csr_matrix`.
    * `gb_zoom_get_stats_for_std_variance_for_sample`: maximum memory required (in GB) for rule `zoom_get_stats_for_std_variance_for_sample`.
    * `gb_zoom_get_mean_var_for_group`: maximum memory required (in GB) for rule `zoom_get_mean_var_for_group`.
    * `gb_zoom_get_estimated_variances`: maximum memory required (in GB) for rule `zoom_get_estimated_variances`.
    * `gb_zoom_get_stats_for_std_variance_for_group`: maximum memory required (in GB) for rule `zoom_get_stats_for_std_variance_for_group`.
    * `gb_zoom_get_highly_variable_genes`: maximum memory required (in GB) for rule `zoom_get_highly_variable_genes`.
    * `gb_zoom_create_hvg_matrix`: maximum memory required (in GB) for rule `zoom_create_hvg_matrix`.
    * `gb_zoom_run_integration`: maximum memory required (in GB) for rule `zoom_run_integration`.
    * `gb_zoom_make_marker_pseudobulks`: maximum memory required (in GB) for rule `zoom_make_marker_pseudobulks`.
    * `gb_zoom_prepare_marker_pseudobulks`: maximum memory required (in GB) for rule `zoom_prepare_marker_pseudobulks`.
    * `gb_zoom_calc_marker_hvgs`: maximum memory required (in GB) for rule `zoom_calc_marker_hvgs`.
    * `gb_zoom_run_marker_genes`: maximum memory required (in GB) for rule `zoom_run_marker_genes`.
    * `gb_zoom_prepare_marker_plot_data`: maximum memory required (in GB) for rule `zoom_prepare_marker_plot_data`.
    * `gb_zoom_run_fgsea`: maximum memory required (in GB) for rule `zoom_run_fgsea`.
    * `gb_zoom_make_subsets`: maximum memory required (in GB) for rule `zoom_make_subsets`.
    * `gb_zoom_train_xgboost_prepare`: maximum memory required (in GB) for rule `zoom_train_xgboost_prepare`.
    * `gb_zoom_train_xgboost_train`: maximum memory required (in GB) for rule `zoom_train_xgboost_train`.
    * `gb_zoom_train_xgboost_predict`: maximum memory required (in GB) for rule `zoom_train_xgboost_predict`.
    * `gb_zoom_train_xgboost_combine`: maximum memory required (in GB) for rule `zoom_train_xgboost_combine`.
    * `gb_zoom_train_xgboost_depth_evaluation`: maximum memory required (in GB) for rule `zoom_train_xgboost_depth_evaluation`.
    * `gb_render_html_zoom_train_xgboost`: maximum memory required (in GB) for rule `render_html_zoom_train_xgboost`.
    * `gb_render_html_zoom`: maximum memory required (in GB) for rule `render_html_zoom`.
    * `mins_run_mapping_flex`: maximum runtime required (in minutes) for rule `run_mapping_flex`.
    * `mins_run_mapping`: maximum runtime required (in minutes) for rule `run_mapping`.
    * `mins_save_alevin_to_h5`: maximum runtime required (in minutes) for rule `save_alevin_to_h5`.
    * `mins_build_hto_index`: maximum runtime required (in minutes) for rule `build_hto_index`.
    * `mins_run_mapping_hto`: maximum runtime required (in minutes) for rule `run_mapping_hto`.
    * `mins_save_alevin_hto_to_h5`: maximum runtime required (in minutes) for rule `save_alevin_hto_to_h5`.
    * `mins_save_alevin_adt_to_h5`: maximum runtime required (in minutes) for rule `save_alevin_adt_to_h5`.
    * `mins_make_hto_sce_objects`: maximum runtime required (in minutes) for rule `make_hto_sce_objects`.
    * `mins_run_cellbender`: maximum runtime required (in minutes) for rule `run_cellbender`.
    * `mins_run_decontx`: maximum runtime required (in minutes) for rule `run_decontx`.
    * `mins_run_cell_calling`: maximum runtime required (in minutes) for rule `run_cell_calling`.
    * `mins_get_barcode_qc_metrics`: maximum runtime required (in minutes) for rule `get_barcode_qc_metrics`.
    * `mins_check_demux_ids`: maximum runtime required (in minutes) for rule `check_demux_ids`.
    * `mins_run_doublet_id_one_run`: maximum runtime required (in minutes) for rule `run_doublet_id_one_run`.
    * `mins_filter_qc_one_run`: maximum runtime required (in minutes) for rule `filter_qc_one_run`.
    * `mins_merge_doublet_diagnostics`: maximum runtime required (in minutes) for rule `merge_doublet_diagnostics`.
    * `mins_merge_qc`: maximum runtime required (in minutes) for rule `merge_qc`.
    * `mins_merge_rowdata`: maximum runtime required (in minutes) for rule `merge_rowdata`.
    * `mins_get_qc_sample_statistics`: maximum runtime required (in minutes) for rule `get_qc_sample_statistics`.
    * `mins_make_one_pb_empty`: maximum runtime required (in minutes) for rule `make_one_pb_empty`.
    * `mins_merge_pb_empty`: maximum runtime required (in minutes) for rule `merge_pb_empty`.
    * `mins_calculate_ambient_genes`: maximum runtime required (in minutes) for rule `calculate_ambient_genes`.
    * `mins_make_one_pb_cells`: maximum runtime required (in minutes) for rule `make_one_pb_cells`.
    * `mins_merge_pb_cells`: maximum runtime required (in minutes) for rule `merge_pb_cells`.
    * `mins_make_tmp_csr_matrix`: maximum runtime required (in minutes) for rule `make_tmp_csr_matrix`.
    * `mins_get_stats_for_std_variance_for_sample`: maximum runtime required (in minutes) for rule `get_stats_for_std_variance_for_sample`.
    * `mins_get_mean_var_for_group`: maximum runtime required (in minutes) for rule `get_mean_var_for_group`.
    * `mins_get_estimated_variances`: maximum runtime required (in minutes) for rule `get_estimated_variances`.
    * `mins_get_stats_for_std_variance_for_group`: maximum runtime required (in minutes) for rule `get_stats_for_std_variance_for_group`.
    * `mins_get_highly_variable_genes`: maximum runtime required (in minutes) for rule `get_highly_variable_genes`.
    * `mins_create_hvg_matrix`: maximum runtime required (in minutes) for rule `create_hvg_matrix`.
    * `mins_create_doublets_hvg_matrix`: maximum runtime required (in minutes) for rule `create_doublets_hvg_matrix`.
    * `mins_calculate_tricycle_projection_one_run`: maximum runtime required (in minutes) for rule `calculate_tricycle_projection_one_run`.
    * `mins_aggregate_tricycle`: maximum runtime required (in minutes) for rule `aggregate_tricycle`.
    * `mins_estimate_tricycle_reference_ellipse`: maximum runtime required (in minutes) for rule `estimate_tricycle_reference_ellipse`.
    * `mins_run_integration`: maximum runtime required (in minutes) for rule `run_integration`.
    * `mins_make_clean_h5ads`: maximum runtime required (in minutes) for rule `make_clean_h5ads`.
    * `mins_convert_h5ad_to_sce`: maximum runtime required (in minutes) for rule `convert_h5ad_to_sce`.
    * `mins_make_marker_pseudobulks`: maximum runtime required (in minutes) for rule `make_marker_pseudobulks`.
    * `mins_prepare_marker_pseudobulks`: maximum runtime required (in minutes) for rule `prepare_marker_pseudobulks`.
    * `mins_calc_marker_hvgs`: maximum runtime required (in minutes) for rule `calc_marker_hvgs`.
    * `mins_run_marker_genes`: maximum runtime required (in minutes) for rule `run_marker_genes`.
    * `mins_prepare_marker_plot_data`: maximum runtime required (in minutes) for rule `prepare_marker_plot_data`.
    * `mins_run_fgsea`: maximum runtime required (in minutes) for rule `run_fgsea`.
    * `mins_render_html_mapping`: maximum runtime required (in minutes) for rule `render_html_mapping`.
    * `mins_render_html_multiplexing`: maximum runtime required (in minutes) for rule `render_html_multiplexing`.
    * `mins_render_html_ambient`: maximum runtime required (in minutes) for rule `render_html_ambient`.
    * `mins_render_html_qc`: maximum runtime required (in minutes) for rule `render_html_qc`.
    * `mins_render_html_hvgs`: maximum runtime required (in minutes) for rule `render_html_hvgs`.
    * `mins_render_html_integration`: maximum runtime required (in minutes) for rule `render_html_integration`.
    * `mins_render_html_marker_genes`: maximum runtime required (in minutes) for rule `render_html_marker_genes`.
    * `mins_render_html_label_celltypes`: maximum runtime required (in minutes) for rule `render_html_label_celltypes`.
    * `mins_render_html_train_xgboost`: maximum runtime required (in minutes) for rule `render_html_train_xgboost`.
    * `mins_run_celltypist`: maximum runtime required (in minutes) for rule `run_celltypist`.
    * `mins_run_scprocess_labeller`: maximum runtime required (in minutes) for rule `run_scprocess_labeller`.
    * `mins_merge_labels`: maximum runtime required (in minutes) for rule `merge_labels`.
    * `mins_train_xgboost_prepare`: maximum runtime required (in minutes) for rule `train_xgboost_prepare`.
    * `mins_train_xgboost_train`: maximum runtime required (in minutes) for rule `train_xgboost_train`.
    * `mins_train_xgboost_predict`: maximum runtime required (in minutes) for rule `train_xgboost_predict`.
    * `mins_train_xgboost_combine`: maximum runtime required (in minutes) for rule `train_xgboost_combine`.
    * `mins_train_xgboost_depth_evaluation`: maximum runtime required (in minutes) for rule `train_xgboost_depth_evaluation`.
    * `mins_build_shiny_data`: maximum runtime required (in minutes) for rule `build_shiny_data`. The pre-split key `mins_build_shiny_app` is also accepted.
    * `mins_build_zoom_shiny_data`: maximum runtime required (in minutes) for rule `build_zoom_shiny_data`. The pre-split key `mins_build_zoom_shiny_app` is also accepted.
    * `mins_zoom_check_clusters`: maximum runtime required (in minutes) for rule `zoom_check_clusters`.
    * `mins_zoom_filter_cells_qc`: maximum runtime required (in minutes) for rule `zoom_filter_cells_qc`.
    * `mins_zoom_aggregate_labels`: maximum runtime required (in minutes) for rule `zoom_aggregate_labels`.
    * `mins_zoom_get_sample_statistics`: maximum runtime required (in minutes) for rule `zoom_get_sample_statistics`.
    * `mins_zoom_make_one_pb_cells`: maximum runtime required (in minutes) for rule `zoom_make_one_pb_cells`.
    * `mins_zoom_merge_pb_cells`: maximum runtime required (in minutes) for rule `zoom_merge_pb_cells`.
    * `mins_zoom_calculate_ambient_genes`: maximum runtime required (in minutes) for rule `zoom_calculate_ambient_genes`.
    * `mins_zoom_make_tmp_csr_matrix`: maximum runtime required (in minutes) for rule `zoom_make_tmp_csr_matrix`.
    * `mins_zoom_get_stats_for_std_variance_for_sample`: maximum runtime required (in minutes) for rule `zoom_get_stats_for_std_variance_for_sample`.
    * `mins_zoom_get_mean_var_for_group`: maximum runtime required (in minutes) for rule `zoom_get_mean_var_for_group`.
    * `mins_zoom_get_estimated_variances`: maximum runtime required (in minutes) for rule `zoom_get_estimated_variances`.
    * `mins_zoom_get_stats_for_std_variance_for_group`: maximum runtime required (in minutes) for rule `zoom_get_stats_for_std_variance_for_group`.
    * `mins_zoom_get_highly_variable_genes`: maximum runtime required (in minutes) for rule `zoom_get_highly_variable_genes`.
    * `mins_zoom_create_hvg_matrix`: maximum runtime required (in minutes) for rule `zoom_create_hvg_matrix`.
    * `mins_zoom_run_integration`: maximum runtime required (in minutes) for rule `zoom_run_integration`.
    * `mins_zoom_make_marker_pseudobulks`: maximum runtime required (in minutes) for rule `zoom_make_marker_pseudobulks`.
    * `mins_zoom_prepare_marker_pseudobulks`: maximum runtime required (in minutes) for rule `zoom_prepare_marker_pseudobulks`.
    * `mins_zoom_calc_marker_hvgs`: maximum runtime required (in minutes) for rule `zoom_calc_marker_hvgs`.
    * `mins_zoom_run_marker_genes`: maximum runtime required (in minutes) for rule `zoom_run_marker_genes`.
    * `mins_zoom_prepare_marker_plot_data`: maximum runtime required (in minutes) for rule `zoom_prepare_marker_plot_data`.
    * `mins_zoom_run_fgsea`: maximum runtime required (in minutes) for rule `zoom_run_fgsea`.
    * `mins_zoom_make_subsets`: maximum runtime required (in minutes) for rule `zoom_make_subsets`.
    * `mins_zoom_train_xgboost_prepare`: maximum runtime required (in minutes) for rule `zoom_train_xgboost_prepare`.
    * `mins_zoom_train_xgboost_train`: maximum runtime required (in minutes) for rule `zoom_train_xgboost_train`.
    * `mins_zoom_train_xgboost_predict`: maximum runtime required (in minutes) for rule `zoom_train_xgboost_predict`.
    * `mins_zoom_train_xgboost_combine`: maximum runtime required (in minutes) for rule `zoom_train_xgboost_combine`.
    * `mins_zoom_train_xgboost_depth_evaluation`: maximum runtime required (in minutes) for rule `zoom_train_xgboost_depth_evaluation`.
    * `mins_render_html_zoom_train_xgboost`: maximum runtime required (in minutes) for rule `render_html_zoom_train_xgboost`.
    * `mins_render_html_zoom`: maximum runtime required (in minutes) for rule `render_html_zoom`.

## {{scjoin}} { #scprocess-join }

**Description**: Integrate multiple completed {{sc}} projects.

**Parameters**:

* `configfile` (positional): path to a configuration YAML file.
* `-n`/`--dry-run` (optional): list the join steps without executing them. Helpful for checking input files and parameters.
* `-E`/`--extraargs` (optional): list of additional arguments to pass to `Snakemake`. Refer to [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli) for a detailed explanation of available command-line options.
* `--unlock` (optional): unlock the directory if a previous run was interrupted.
* `-r`/`--rule` (optional): specifies which rule {{sc}} should run. Options are:
    + `all`: default.
    + `label_celltypes`: cell type annotation.
    + `shiny`: interactive Shiny app generation.
    + `train_xgboost`: `XGBoost` classifier training.
    + `zoom`: subclustering.
* `--zoom` (optional): with `-r zoom` or `-r shiny`, select a named joined zoom, or `all`. Without this option, `-r zoom` runs every listed zoom; `-r shiny` builds the parent join's app.

### configuration file

The config file for {{scjoin}} requires `join` and `projects` section and can optionally contain `hvg`, `integration`, `marker_genes`, `label_celltypes`, `train_xgboost`, and `shiny` sections with same parameters as the {{scrun}} config file, with the following exceptions:

* `hvg`: `hvg_method`, `hvg_chunk_size`, `hvg_exclude_ambient_genes`, `hvg_exclude_from_file`,and `hvg_metadata_split_var` are not available
* `integration`: `int_dbl_res`, `int_dbl_cl_prop`, `int_dbl_threshold_method` and `int_sce_outs` are not available
* `label_celltypes`: `min_cl_size` is not available
* `marker_genes`: `mkr_adjust_project_id` is specific to {{scjoin}}: boolean specifying whether the marker gene model fits a `project_id` term alongside `cluster`, correcting for systematic baseline differences between source projects when computing markers from the joint pseudobulks. Default: `false`.

This is an example config file for {{scjoin}}. Required parameters are highlighted:

```yaml hl_lines="1 2 3 4 5 6 7 9"
join:
  name: my_join
  proj_dir: /path/to/join_output
  date_stamp: "2025-01-15"
  ref_txome: human_2024 # For Flex data, replace ref_txome with probe_set:
  your_name: Testy McUser
  affiliation: where you work
  metadata_vars: [var1, var2, var3]
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
# arrays can be used for `int_batch_var` and `int_theta` to configure multiple Harmony covariates. When both values are arrays, they must have the same length so that each batch variable has a corresponding theta value.
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

#### join

* `name`: name of the joint project.
* `proj_dir`: absolute path to `workflowr` project directory created with the {{scnewjoin}} function.
* `your_name`: author’s name, displayed in HTML outputs.
* `affiliation`: author’s affiliation, displayed in HTML outputs.
* `date_stamp`: start date of the analysis, formatted as `"YYYY-MM-DD"`.
* `ref_txome`: referece transcriptome used for mapping source projects; Required for joining polyA data; only compatible projects can be joined.
* `probe_set`: probe set used for mapping source projects; Required for joining Flex data; only compatibe projects can be joined.
* `metadata_vars` (optional): list of column names from the source projects' `sample_metadata` CSV files to carry through into the joint outputs. Each variable must exist in at least one project's `sample_metadata` file.

#### projects

Mapping containing at least two source projects. Each entry requires `config` - the path to a completed project config, and may specify `zoom_name` to join a completed zoom instead of the full project.

#### resources  

This section allows users to adjust the resource requirements for specific `Snakemake` rules in {{scjoin}}. By specifying the parameters below, users can fine-tune these settings for their pipeline:

* `gb_[rule_name]`: specifies the maximum memory (in GB) requested for running a specific rule. `rule_name` should be replaced with an {{sc}} rule name. This value applies for the entire job, not per thread.
* `mins_[rule_name]`: specifies the maximum runtime (in minutes) requested for running a specific rule. `rule_name` should be replace with an {{sc}} rule name.

Additional parameters include:

* `retries`: number of times to retry running a specific rule in {{scjoin}} if it fails. For each attempt the initial memory requested for the rule is multiplied by `1.5**(attempt - 1)`. Useful for when {{scjoin}} is ran on a [cluster](setup.md#cluster-setup).
* `n_join_marker_genes_calc`: number of threads requested for rule `join_marker_genes_calc`. Default is 8.

??? note "Detailed information about resource parameters"

    * `gb_join_select_hvgs`: maximum memory required (in GB) for rule `join_select_hvgs`.
    * `gb_join_build_matrix`: maximum memory required (in GB) for rule `join_build_matrix`.
    * `gb_join_build_coldata`: maximum memory required (in GB) for rule `join_build_coldata`.
    * `gb_join_build_sample_metadata`: maximum memory required (in GB) for rule `join_build_sample_metadata`.
    * `gb_join_build_qc`: maximum memory required (in GB) for rule `join_build_qc`.
    * `gb_join_build_bpcells_matrix`: maximum memory required (in GB) for rule `join_build_bpcells_matrix`.
    * `gb_join_pca`: maximum memory required (in GB) for rule `join_pca`.
    * `gb_join_integration`: maximum memory required (in GB) for rule `join_integration`.
    * `gb_join_build_h5ads_yaml`: maximum memory required (in GB) for rule `join_build_h5ads_yaml`.
    * `gb_join_make_pseudobulks`: maximum memory required (in GB) for rule `join_make_pseudobulks`.
    * `gb_join_marker_genes_prep`: maximum memory required (in GB) for rule `join_marker_genes_prep`.
    * `gb_join_calc_hvgs`: maximum memory required (in GB) for rule `join_calc_hvgs`.
    * `gb_join_marker_genes_calc`: maximum memory required (in GB) for rule `join_marker_genes_calc`.
    * `gb_join_marker_plot_data`: maximum memory required (in GB) for rule `join_marker_plot_data`.
    * `gb_join_fgsea`: maximum memory required (in GB) for rule `join_fgsea`.
    * `gb_join_celltypist`: maximum memory required (in GB) for rule `join_celltypist`.
    * `gb_join_scprocess_labeller`: maximum memory required (in GB) for rule `join_scprocess_labeller`.
    * `gb_join_extract_labels`: maximum memory required (in GB) for rule `join_extract_labels`.
    * `gb_join_merge_labels`: maximum memory required (in GB) for rule `join_merge_labels`.
    * `gb_join_render_html`: maximum memory required (in GB) for rule `join_render_html`.
    * `gb_join_train_xgboost_prepare`: maximum memory required (in GB) for rule `join_train_xgboost_prepare`.
    * `gb_join_train_xgboost_train`: maximum memory required (in GB) for rule `join_train_xgboost_train`.
    * `gb_join_train_xgboost_predict`: maximum memory required (in GB) for rule `join_train_xgboost_predict`.
    * `gb_join_train_xgboost_combine`: maximum memory required (in GB) for rule `join_train_xgboost_combine`.
    * `gb_join_train_xgboost_depth_evaluation`: maximum memory required (in GB) for rule `join_train_xgboost_depth_evaluation`.
    * `gb_join_render_html_train_xgboost`: maximum memory required (in GB) for rule `join_render_html_train_xgboost`.
    * `gb_build_shiny_data`: maximum memory required (in GB) for rule `build_shiny_data`. The pre-split key `gb_build_shiny_app` is also accepted.
    * `gb_build_zoom_shiny_data`: maximum memory required (in GB) for rule `build_zoom_shiny_data`. The pre-split key `gb_build_zoom_shiny_app` is also accepted.
    * `gb_join_zoom_prepare`: maximum memory required (in GB) for rule `join_zoom_prepare`.
    * `gb_join_zoom_calculate_raw_gene_count_statistics`: maximum memory required (in GB) for rule `join_zoom_calculate_raw_gene_count_statistics`.
    * `gb_join_zoom_fit_expected_gene_variances`: maximum memory required (in GB) for rule `join_zoom_fit_expected_gene_variances`.
    * `gb_join_zoom_calculate_clipped_gene_count_statistics`: maximum memory required (in GB) for rule `join_zoom_calculate_clipped_gene_count_statistics`.
    * `gb_join_zoom_rank_and_select_highly_variable_genes`: maximum memory required (in GB) for rule `join_zoom_rank_and_select_highly_variable_genes`.
    * `mins_join_select_hvgs`: maximum runtime required (in minutes) for rule `join_select_hvgs`.
    * `mins_join_build_matrix`: maximum runtime required (in minutes) for rule `join_build_matrix`.
    * `mins_join_build_coldata`: maximum runtime required (in minutes) for rule `join_build_coldata`.
    * `mins_join_build_sample_metadata`: maximum runtime required (in minutes) for rule `join_build_sample_metadata`.
    * `mins_join_build_qc`: maximum runtime required (in minutes) for rule `join_build_qc`.
    * `mins_join_build_bpcells_matrix`: maximum runtime required (in minutes) for rule `join_build_bpcells_matrix`.
    * `mins_join_pca`: maximum runtime required (in minutes) for rule `join_pca`.
    * `mins_join_integration`: maximum runtime required (in minutes) for rule `join_integration`.
    * `mins_join_build_h5ads_yaml`: maximum runtime required (in minutes) for rule `join_build_h5ads_yaml`.
    * `mins_join_make_pseudobulks`: maximum runtime required (in minutes) for rule `join_make_pseudobulks`.
    * `mins_join_marker_genes_prep`: maximum runtime required (in minutes) for rule `join_marker_genes_prep`.
    * `mins_join_calc_hvgs`: maximum runtime required (in minutes) for rule `join_calc_hvgs`.
    * `mins_join_marker_genes_calc`: maximum runtime required (in minutes) for rule `join_marker_genes_calc`.
    * `mins_join_marker_plot_data`: maximum runtime required (in minutes) for rule `join_marker_plot_data`.
    * `mins_join_fgsea`: maximum runtime required (in minutes) for rule `join_fgsea`.
    * `mins_join_celltypist`: maximum runtime required (in minutes) for rule `join_celltypist`.
    * `mins_join_scprocess_labeller`: maximum runtime required (in minutes) for rule `join_scprocess_labeller`.
    * `mins_join_extract_labels`: maximum runtime required (in minutes) for rule `join_extract_labels`.
    * `mins_join_merge_labels`: maximum runtime required (in minutes) for rule `join_merge_labels`.
    * `mins_join_render_html`: maximum runtime required (in minutes) for rule `join_render_html`.
    * `mins_join_train_xgboost_prepare`: maximum runtime required (in minutes) for rule `join_train_xgboost_prepare`.
    * `mins_join_train_xgboost_train`: maximum runtime required (in minutes) for rule `join_train_xgboost_train`.
    * `mins_join_train_xgboost_predict`: maximum runtime required (in minutes) for rule `join_train_xgboost_predict`.
    * `mins_join_train_xgboost_combine`: maximum runtime required (in minutes) for rule `join_train_xgboost_combine`.
    * `mins_join_train_xgboost_depth_evaluation`: maximum runtime required (in minutes) for rule `join_train_xgboost_depth_evaluation`.
    * `mins_join_render_html_train_xgboost`: maximum runtime required (in minutes) for rule `join_render_html_train_xgboost`.
    * `mins_build_shiny_data`: maximum runtime required (in minutes) for rule `build_shiny_data`. The pre-split key `mins_build_shiny_app` is also accepted.
    * `mins_build_zoom_shiny_data`: maximum runtime required (in minutes) for rule `build_zoom_shiny_data`. The pre-split key `mins_build_zoom_shiny_app` is also accepted.
    * `mins_join_zoom_prepare`: maximum runtime required (in minutes) for rule `join_zoom_prepare`.
    * `mins_join_zoom_calculate_raw_gene_count_statistics`: maximum runtime required (in minutes) for rule `join_zoom_calculate_raw_gene_count_statistics`.
    * `mins_join_zoom_fit_expected_gene_variances`: maximum runtime required (in minutes) for rule `join_zoom_fit_expected_gene_variances`.
    * `mins_join_zoom_calculate_clipped_gene_count_statistics`: maximum runtime required (in minutes) for rule `join_zoom_calculate_clipped_gene_count_statistics`.
    * `mins_join_zoom_rank_and_select_highly_variable_genes`: maximum runtime required (in minutes) for rule `join_zoom_rank_and_select_highly_variable_genes`.



