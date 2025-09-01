# Reference

## {{scsetup}} { #scprocess-setup }

**Description**: Download all data required for {{sc}} and index reference genomes for `simpleaf`.

**Parameters**:
The command requires a configuration file named `scprocess_setup.yaml` located in {{sc}} data directory (for instructions on how to set up the {{sc}} data directory see the [Getting started](setup.md#scprocess-data-directory-setup) section). In this file, the user has to specify which reference genome will be made available for {{sc}}. For example:

```yaml
genome:
  tenx:
    - name: human_2024 
      decoys: True
      rrnas: True
  custom:
    - name: custom_genome_name
      fasta: /path/to/genome.fa
      gtf: /path/to/genes.gtf
      decoys: True
      mito_str: "^mt-"
    - name: custom_genome_name2
      index_dir: /path/to/prebuild/alevin/index
      gtf: /path/to/genes.gtf
      mito_str: "^MT-"
```

Prebuilt human and mouse reference genomes from 10x Genomics can be downloaded with {{scsetup}} by adding `tenx` to the `scprocess_setup.yaml` file. Valid values for names are `human_2024`, `mouse_2024`, `human_2020`, `mouse_2020`.  

Names and specifications for custom references should be listed in the `custom` section of the `scprocess_setup.yaml` file. For each `custom` genome users have to provide the following parameters:

* one of:
    + `fasta`: path to FASTA file
    + `index_dir`: path to prebuild alevin index; when specified `decoys` option is ignored
* `gtf`: path to GTF file [(specific format?)]
* `mito_str`: regular expression used to identify genes in the mitochondial genome (example for mouse: `"^mt-"`)

Optional parameters for both `tenx` and `custom` references are:

* `decoys`: whether or not poison k-mer information should be inserted into the index. This parameter is optional. If not specified, it defaults to `True` for all genomes.

Optional paramater for `tenx` references is:

* `rrnas`: whether or not ribosomal RNAs should be included in the reference. If not specified it defaults to `True` for all `tenx` genomes.

!!! note "Impact of custom parameters for `tenx` genomes on `scsetup` runtime"

    When configuring `tenx` genomes with their default values, `scsetup` will download prebuilt indices optimized for `simpleaf`. However, if the default parameters are modified (e.g., setting `rrnas` or `decoys` to `False`), `scsetup` will build the indices from scratch during execution, which will increase the runtime.


!!! info "More about decoys"
    {{sc}} utilizes `simpleaf`, a lightweight mapping approach that, by default, maps sequenced fragments exclusively to the transcriptome. However, this can lead to incorrect mapping of reads that arise from unannotated genomic loci to the transcriptome. To mitigate this issue, the `decoys` parameter in `scsetup` is set to `True`. This option allows `simpleaf` to identify genomic regions with sequences similar to those in transcribed regions (decoys), thereby reducing the likelihood of false mappings. We strongly recommend keeping the decoy setting enabled. For further details, refer to [Srivastava et al., 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8).


## {{scnew}}

**Description**: Create a new `workflowr` project directory for {{sc}} outputs.

**Parameters**: 

* `name` (positional): name of the new `workflowr` project directory
* `-w`/`--where` (optional): path to the directory where the new project will be created; defaults to the current working directory
* `-s`/`--sub` (optional): if provided, creates `data/fastqs` and `data/metadata` subdirectories within the project.
* `-c`/`--config` (optional): if provided, generates a template configuration YAML file for {{sc}}


## {{scrun}}

**Description**: Run {{sc}}.

**Parameters**:

* `-n`/`--dry-run`: This makes {{sc}} perform a trial run that doesn't create any new files. This is helpful for (1) checking that the various input files and parameter settings are likely to work, and (2) checking what work will be done by {{sc}}.
* `-E`/`--extraagrs`: `snakemake` is a sophisticated package with many options that can be set by the user. This argument allows users to set additional arguments; see [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for `snakemake`'s documentation of command line arguments.
* `-r`/`--rule`: Specifies which rule {{sc}} should run. The options are:
    + `all`: default; includes all [Core pipeline steps](introduction.md#core-pipeline-steps)
    + `mapping`: read alignment and quantification using `simpleaf`.
    + `cellbender`: ambient RNA removal with cellbender.
    + `demux`: sample demultiplexing.
    + `qc`: qc filtering.
    + `hvg`: calculation of highly variable genes.
    + `integration`: dimentionality reduction with PCA and optional batch correction with `Harmony`.
    + `marker_genes`: marker gene identification.
    + `label_celltypes`: cell type annotation using a pre-trained classifier or a custom annotation file. 


## {{scknee}} { #scprocess-plotknee }

**Description**: Create an interactive barcode-rank plot. Can only be used once the mapping step is completed.

**Parameters**:

Specify either `-k`/`--kneefile` or `-c`/`--configfile`:

* `-k`/`--kneefile`: path to the knee plot data file generated by {{sc}}, e.g. `output/[short_tag]_alevin_fry/af_[sample_id]/knee_plot_data_[sample_id]_[date_stamp].txt.gz`.
* `-c`/`--configfile`: path to configuration file used for running {{sc}}
* `-s`/`--sample`: sample_id corresponding to barcode-rank curve


### configuration file

This is an example config file for {{sc}} with all parameters and their default values/placeholders. Required parameters are highlighted:

=== "default values"

    ```yaml hl_lines="1 2 3 4 5 6 7 8 9 10 11"
    proj_dir: 
    fastq_dir:  
    full_tag: 
    short_tag:
    your_name:
    affiliation:
    sample_metadata:
    species:
    date_stamp:
    alevin:
      chemistry:
    custom_sample_params:
    exclude:
      sample_id:
      pooll_id:
    metadata_vars:
    multiplexing:
    ambient:
      ambient_method: decontx
      cellbender_version: 'v0.3.0'
      cell_calling: barcodeRanks
      cb_max_prop_kept: 0.9 
      cb_force_expected_cells:
      cb_force_total_droplets_included:
      cb_force_low_count_threshold:
      cb_force_learning_rate:
    qc:
      qc_min_counts: 300   
      qc_min_feats: 300        
      qc_min_mito: 0             
      qc_max_mito: 0.1         
      qc_min_splice: 0          
      qc_max_splice: 1          
      qc_min_cells: 500 
      dbl_min_feats: 100
    hvg:
      hvg_method: sample
      hvg_metadata_split_var: 
      n_hvgs: 2000
      exclude_empties: True
    integration:  
      cl_method: louvain
      reduction: harmony
      int_n_dims:       50                 
      int_dbl_res:      4                   
      int_dbl_cl_prop:  0.5                 
      int_theta:        0.1                 
      int_res_ls:       [0.1, 0.2, 0.5, 1, 2]  
    marker_genes:
      custom_sets:
        - name: 
          file: 
      mkr_min_cl_size: 100                                  
      mkr_min_cells: 10                                      
      mkr_not_ok_re: "(lincRNA|lncRNA|pseudogene|antisense)"  
      mkr_min_cpm_mkr: 50                                     
      mkr_min_cpm_go: 1                                      
      mkr_max_zero_p: 0.5                                    
      mkr_gsea_cut: 0.1
      mkr_sel_res: 0.2
    label_celltypes:
      lbl_tissue:      
      lbl_sel_res_cl:  "RNN_snn_res.2"  
      lbl_min_pred:    0.5              
      lbl_min_cl_prop: 0.5             
      lbl_min_cl_size: 100
    zoom:
    resources:
      mb_run_mapping: 8192               
      mb_save_alevin_to_h5: 8192            
      mb_run_ambient: 32768            
      mb_run_hvgs: 8192  
      mb_run_scdblfinder: 4096             
      mb_combine_scdblfinder_outputs: 8192  
      mb_run_harmony: 16384                 
      mb_run_marker_genes: 16384            
      mb_html_marker_genes: 8192            
      mb_lbl_label_celltypes: 16384         
      mb_lbl_save_subset_sces: 16384        
      mb_lbl_render_template_rmd: 4096
    ```

=== "placeholders"

    ```yaml hl_lines="1 2 3 4 5 6 7 8 9 10 11"
    proj_dir: /path/to/proj/directory #
    fastq_dir: /path/to/directory/with/fastq/files 
    full_tag: test_project 
    short_tag: test 
    your_name: John Doe 
    affiliation: where you work 
    sample_metadata: /path/to/metadata.csv
    species: human_2024
    date_stamp: "2050-01-01" 
    alevin:
      chemistry: 3v3
    custom_sample_params: /path/to/file/with/custom_parameters.yaml # modify
    exclude:
      sample_id: [sample1, sample2]
      pool_id: [pool1, pool2]
    metadata_vars: [test1, test2]
    ambient:
      ambient_method: decontx
      cellbender_version: 'v0.3.0'
      cell_calling: barcodeRanks
      cb_max_prop_kept: 0.9 
      cb_force_expected_cells: 10000
      cb_force_total_droplets_included: 20000
      cb_force_low_count_threshold: 5
      cb_force_learning_rate: 0.001
    multiplexing:
      demux_type: af
      fastq_dir: /path/to/directory/with/hto_fastq/files 
      feature_ref: /path/to/feature_ref.csv
      batch_var: sample_id
    qc:
      qc_min_counts: 300   
      qc_min_feats: 300        
      qc_min_mito: 0             
      qc_max_mito: 0.1         
      qc_min_splice: 0          
      qc_max_splice: 1          
      qc_min_cells: 500  
      dbl_min_feats: 100
    hvg:
      hvg_method: sample
      n_hvgs: 2000
      exclude_empties: True
    integration:
      cl_method: louvain
      reduction: harmony    
      int_n_dims:       50                 
      int_dbl_res:      4                   
      int_dbl_cl_prop:  0.5                 
      int_theta:        0.1                 
      int_res_ls:       [0.1, 0.2, 0.5, 1, 2]
    marker_genes:
      custom_sets:
        - name: mouse_brain
          file: /path/to/file/with/marker/genes.csv
      mkr_min_cl_size: 100                                  
      mkr_min_cells: 10                                      
      mkr_not_ok_re: "(lincRNA|lncRNA|pseudogene|antisense)"  
      mkr_min_cpm_mkr: 50                                     
      mkr_min_cpm_go: 1                                      
      mkr_max_zero_p: 0.5                                    
      mkr_gsea_cut: 0.1
      mkr_sel_res: 0.2
    label_celltypes:
      lbl_tissue:      "brain_cns"      
      lbl_sel_res_cl:  "RNN_snn_res.2"  
      lbl_min_pred:    0.5              
      lbl_min_cl_prop: 0.5             
      lbl_min_cl_size: 100              
    zoom:
      cell_subset_1: /path/to/cell_subset_1_zoom_params.yaml
      cell_subset_2: /path/to/cell_subset_2_zoom_params.yaml
      cell_subset_3: /path/to/cell_subset_3_zoom_params.yaml
    resources:
      mb_run_mapping: 8192               
      mb_save_alevin_to_h5: 8192            
      mb_run_ambient: 32768            
      mb_run_hvgs: 8192 
      mb_run_scdblfinder: 4096             
      mb_combine_scdblfinder_outputs: 8192            
      mb_run_harmony: 16384                 
      mb_run_marker_genes: 16384            
      mb_html_marker_genes: 8192            
      mb_lbl_label_celltypes: 16384         
      mb_lbl_save_subset_sces: 16384        
      mb_lbl_render_template_rmd: 4096
    ```


#### Required parameters

* `proj_dir`: absolute path to `workflowr` project directory created with the `newproj` function.
* `fastq_dir`: path to directory containing FASTQ files. Should be absolute or relative to `proj_dir`.
* `full_tag`: full project label, used in output file names.
* `short_tag`: abbreviated project label, used in output directory names.
* `your_name`: author’s name, displayed in HTML outputs.
* `affiliation`: author’s affiliation, displayed in HTML outputs.
* `sample_metadata`: path to CSV file with sample metadata. Should be absolute or relative to `proj_dir`. Spaces in column names are not allowed. Only required column is `sample_id`; values in `sample_id` should not contain `_R_`and `_R2_`strings.
* `species`: must match one of the values in the `genome_name` column of `setup_parameters.csv` (created by `scsetup`).
* `date_stamp`: start date of the analysis, formatted as `"YYYY-MM-DD"`.
* `chemistry`: 10x assay configurtaion. Accepted values are `3LT`, `3v2`, `3v3`, `3v4`, `5v1`, `5v2`, `5v3`, and `multiome`. `multiome` refers only to gene expression data genertaed with the 10x multiome kit (ATACseq data is not supported).

#### Optional parameters

##### general
 
* `custom_sample_params`: YAML file with optional custom parameters for each sample (custom chemistry, custom ambient and custom cellbender parameters can be specified for each sample). Example:

```yaml
sample_1:
  chemistry: 5v2
  ambient:
    knee1: 4000
    shin1: 400
    knee2: 30
    shin2: 5
sample_2: 
  chemistry: 5v2
  ambient:
    knee1: 3000
    shin1: 400
    knee2: 30
    shin2: 5
sample_3: 
  cellbender:
    total_droplets_included: 20000
```

* `exclude`: list of all samples that should be excluded from the analysis. Samples can be listed under `pool_id`(if multiplexed) or `sample_id`
* `metadata_vars`: list of all variables in `sample_metadata` file. Enables plotting of cell proportions associated with specific metadata values per cluster and in binned UMAPs, with facets representing each variable’s unique values to help assess whether cells with certain annotations are concentrated in specific clusters/UMAP regions, or if they are evenly distributed.


##### multiplexing

* `demux_type`: `af` if demultiplexing of samples should be performed with {{sc}} or `custom` if demultiplexing results will be used as input to {{sc}}
* `fastq_dir`: path to directory containing HTO FASTQ files. Should be absolute or relative to `proj_dir`. Required if `demux_type` is `af`. 
* `feature_ref`: path to CSV file with columns `hto_id` and `sequence`. Required if `demux_type` is `af`.
* `demux_output`: path to CSV file with columns `pool_id`, `sample_id`, `cell_id`. Optional column `class` can be added with values `doublet`, `singlet` or `ambiguous`. Required if `demux_type` is `af`. 
* `batch_var`: variable to use for integration with `harmony`. Options are `pool_id` or `sample_id`. 

##### ambient

* `ambient_method`: method for ambient RNA removal; options are `decontx` (default), `cellbender` or `none`.
* `cellbender_version`: version of `cellbender` to use if `ambient_method` is set to `cellbender`. Options are `'v0.3.0'` (default) and `'v0.2.0'`.
* `cb_force_expected_cells`: forces the `--expected-cells` parameter to be consistent across all samples; applicable only if `ambient_method` is `cellbender`. For more information about this parameter see [Cellbender's documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_force_total_droplets_included`: forces the `--total-droplets-included` parameter to be consistent across all samples; applicable only if `ambient_method` is `cellbender`. For more information about this parameter see [Cellbender's documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_force_low_count_threshold`: forces the `--low-count-threshold` parameter to be consistent across all samples; applicable only if `ambient_method` is `cellbender`. For more information about this parameter see [Cellbender's documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_force_learning_rate`: Forces the `--learning-rate` parameter to be consistent across all samples; applicable only if `ambient_method` is `cellbender`. For more information about this parameter see [Cellbender's documentation](https://cellbender.readthedocs.io/en/latest/reference/index.html).
* `cb_max_prop_kept`: maximum proportion of droplets, relative to `--total-droplets-included`, that `cellbender` can call as cells. Default is `0.9`, meaning samples are excluded if `cellbender` calls more than 90% of `--total-droplets-included` droplets as cells. Applicable only if `ambient_method` is `cellbender`.
* `cell_calling`: method for cell calling when `ambient_method` is `none` or `decontx`. Options are `barcodeRanks` and `emptyDrops`.

##### qc

* `qc_min_counts`: minimum number of UMIs per cell required to retain the cell.
* `qc_min_feats`: minimum number of features per cell required to retain the cell.
* `qc_min_mito`: minimum proportion of mitochondrial reads required to retain the cell.
* `qc_max_mito`: maximum proportion of mitochondrial reads allowed to retain the cell.
* `qc_min_splice`: minimum proportion of spliced reads required to retain the cell.
* `qc_max_splice`: maximum proportion of spliced reads allowed to retain the cell.
* `qc_min_cells`: minimum number of cells required in a sample after QC filtering to retain the sample.
* `dbl_min_feats`: number of features required for each barcode to be included in `scDblFinder` calculations.
* `exclude_mito`: boolean; whether to exclude mitochondrial genes or not

##### pb_empties

* `ambient_genes_logfc_thr`: logfc treshold for filtering edgeR results, default: 0
* `ambient_genes_fdr_thr`: fdr threshold for filtering edgeR results, default: 0.01

##### hvg

* `hvg_method`: options: 
    + `sample` - calculate highly variable genes per sample, then calculate combined ranking across samples; 
    + `all` - calculate highly variable genes across all cells in the dataset; and
    + `group` - calculate highly variable genes for each sample group then calculate combined ranking across groups.
* `hvg_metadata_split_var`: if `hvg_method` is `group`, which variable in `sample_metadata` should be used to define sample groups.
* `n_hvgs`: number of HVGs to use for PCA
* `hvg_exclude_ambient_genes`: if `True`, genes enriched in "empty" droplets relative to cells will be excluded from highly variable genes selection. 

##### integration

* `cl_method`: algorithm used for clustering, options: `leiden`, `louvain`
* `reduction`: which dimensional reduction to use for clustering and UMAP, options: `pca`, `harmony`
* `int_n_dims`: number of principal components to use for data integration.
* `int_dbl_res`: clustering resolution for identification of additional doublets.
* `int_dbl_cl_prop`: proportion threshold of doublets in a cluster; clusters exceeding this proportion are excluded.
* `int_theta`: theta parameter for `Harmony` integration, controlling batch variable mixing. `0` means no extra mixing of batch variable; Default in {{sc}} is `0.1`, otherwise `2`. 
* `int_res_ls`: list of cluster resolutions for `Harmony`-based clustering.

##### marker_genes

* `custom_sets`: a list of custom marker gene sets, each defined by a unique name and associated file path.
    + `name`: a string representing the name of the marker gene set
    + `file`: path to CSV file containing a list of genes in the marker gene set. Must contain column `label` (marker gene category), and `symbol` and/or `ensembl_id`. If not speficied `scprocess` will look for file `$SCPROCESS_DATA_DIR/marker_genes/{name}.csv`
* `mkr_min_cl_size`: minimum number of cells required in a cluster to calculate marker genes for that cluster.
* `mkr_min_cells`: minimum number of cells required in a pseudobulk sample to include it in marker gene calculations.
* `mkr_not_ok_re`: regular expression pattern to exclude specific gene types from plots showing marker gene expression. 
* `mkr_min_cpm_mkr`: minimum counts per million (CPM) in a cell type required for a gene to be considered a marker gene.
* `mkr_min_cpm_go`: minimum counts per million (CPM) in a cell type required for a gene to be used in Gene Ontology (GO) analysis.
* `mkr_max_zero_p`: maximum proportion of pseudobulk samples for a cell type that can have zero counts for a gene to be used in GO analysis.
* `mkr_gsea_cut`: False discovery rate (FDR) cutoff for Gene Set Enrichment Analysis (GSEA).
* `mkr_sel_res`: selected cluster resolution used for identifying marker genes.

##### label_celltypes

* `lbl_tissue`: target tissue for cell type labeling. Options are `human_cns`, `mouse_cns`, `human_pbmc`, and `mouse_pbmc` (Only `human_cns` is available at the moment)
* `lbl_sel_res_cl`: selected cluster resolution for cell type labeling (must match one of the values in `int_res_ls`); higher values are recommended for optimal performance.
* `lbl_min_pred`: minimum probability threshold for assigning a cell to a cell type.
* `lbl_min_cl_prop`: minimum proportion of cells in a cluster that need to be labeled for that cluster to be labeled.
* `lbl_min_cl_size`: minimum number of cells in a cluster required for that cluster to be labeled.

##### zoom

Name of each cell subset should be specified, followed by the path to a corresponding YAML file containing the parameters for that subset. Some parameters in the YAML file inherit their definitions from the primary {{sc}} configuration file, including `hvg_method`, `hvg_metadata_split_var`, `n_hvgs`, `hvg_exclude_ambient_genes`, `ambient_genes_logfc_thr`, `ambient_genes_fdr_thr`, `cl_method`, `reduction`, `int_n_dims`, `int_theta`, `int_res_ls`, `int_sel_res`, `custom_sets`, `mkr_min_cl_size`, `mkr_min_cells`, `mkr_not_ok_re`, `mkr_min_cpm_mkr`, `mkr_min_cpm_go`, `mkr_max_zero_p`, `mkr_gsea_cut`, and `mkr_sel_res`

Additional parameters include:

* `labels`: a list of all labels that define cell types/clusters to be included in subclustering (required).
* `labels_source`: specifies how a cell subset is defined (required). Options include:
    - `xgboost`: labels assigned by the XGBoost classifier (using rule `label_celltypes`)
    - `clusters`: labels based on clustering results obtained with {{sc}}
    - `custom`: user-defined cell type annotations
* `min_n_sample`: the minimum number of cells that a sample must have, after subsetting, to remain in the analysis.
* `make_subset_sces`: whether to create SingleCellExperiment objects containing cells that have been assigned one of the values in `labels`.
* `cluster_res`: required if `labels_source` is set to `clusters`; selected clustering resolution values matching one of the values in `int_sel_res`.
* `custom_labels_f`: required if `labels_source` is set to custom; path to CSV file with columns `sample_id`, `cell_id` and `label`.
* `lbl_sel_res_cl`: equivalent of `lbl_sel_res_cl` parameter used for rule `label_celltypes`. Applicable only if `labels_source` is set to `xgboost`.

 
##### resources

* `retries`: number of times to retry running a specific rule in {{sc}} if it fails. For each attempt, the memory requirement for the rule increases by multiplying the base memory by the attempt number. Useful for when {{sc}} is used on a [cluster](setup.md#cluster-setup).
* `mb_run_mapping`: maximum memory required (in MB) for running `simpleaf`. Value applies to the entire job, not per thread.
* `mb_save_alevin_to_h5`:  maximum memory required (in MB) to save `simpleaf` output to H5 format. Value applies to the entire job, not per thread.
* `mb_run_ambient`: maximum memory required (in MB) to run the ambient RNA removal step. Value applies to the entire job, not per thread.
* `mb_get_barcode_qc_metrics`: maximum memory required (in MB) to obtain quality control metrics related to ambient RNA removal. Value applies to the entire job, not per thread.
* `mb_run_scdblfinder`: maximum memory required (in MB) to run `scDblFinder` for doublet detection. Value applies to the entire job, not per thread.
* `mb_combine_scdblfinder_outputs`: maximum memory required (in MB) to combine `scDblFinder` outputs across samples. Value applies to the entire job, not per thread.
* `mb_make_sce_object`: maximum memory required (in MB) to create a `SingleCellExperiment` object. Value applies to the entire job, not per thread.
* `mb_run_harmony`: maximum memory required (in MB) to run `Harmony` integration for batch correction. Value applies to the entire job, not per thread.
* `mb_run_marker_genes`: maximum memory required (in MB) for marker gene identification. Value applies to the entire job, not per thread.
* `mb_lbl_label_celltypes`: maximum memory required (in MB) to label cell types. Value applies to the entire job, not per thread.
* `mb_lbl_save_subset_sces`: maximum memory required (in MB) to save `SingleCellExperiment` objects with cell subsets. Value applies to the entire job, not per thread.
* `mb_lbl_render_template_rmd`: maximum memory required (in MB) to render HTML document from markdown template. Value applies to the entire job, not per thread.
