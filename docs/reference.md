# Reference

## Functions

### `scsetup`

### `newproj`

### {{ software_name }}

### `plotKnee`

## Configuration files

### .scprocess_setup.yaml

```yaml
genome:
  tenx:
    - name: human_2024 
  custom:
    - name: mouse
      fasta: '/path/to/genome.fa'
      gtf: '/path/to/genes.gtf'
      decoys: False
      mito_str: "^mt-"
```

tenx: valid values are `human_2024`, `mouse_2024`, `human_2020` and `mouse_2020` 

#### 10x Genomics prebuild references
    
Human and mouse reference genomes from 10x Genomics can be downloaded with `scsetup` by adding `tenx` to the `.scprocess_setup.yaml` file. Valid values for names are `human_2024`, `mouse_2024`, `human_2020`, `mouse_2020`.  

#### Custom references

* `fasta`: path to FASTA file 
* `gtf`: path to GTF file that has to contain columns `Feature`, `gene_id`, `gene_name`, `gene_type`, `Chromosome`, `Start`, `End`, `Strand` [mabye something else? is the format of these files always the same?]
* `decoys`: whether or not poison k-mer information should be inserted into the index. This parameter is optional. If not specified, it defaults to `True` for all genomes. 
* `mito_str`: regular expression used to identify genes in the mitochondial genome (example for mouse: "^mt-")

!!! info "More about decoys"
    {{ software_name }} utilizes simpleaf, a lightweight mapping approach that, by default, maps sequenced fragments exclusively to the transcriptome. However, this can lead to incorrect mapping of reads that arise from unannotated genomic loci to the transcriptome. To mitigate this issue, the `decoys` parameter in `scsetup` it to `True`. This option allows simpleaf to identify genomic regions with sequences similar to those in transcribed regions (decoys), thereby reducing the likelihood of false mappings. We strongly recommend keeping the decoy setting enabled. For further details, refer to [Srivastava et al., 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8).


### scprocess config file

This is an example config file for {{ software_name }} with all parameters and their default values:

```yaml

proj_dir: /path/to/proj/directory
fastq_dir: /path/to/directory/with/fastq/files
full_tag: test_project_204
short_tag: test
your_name: John Doe
affiliation: where you work
sample_metadata: /path/to/metadata
species: human_2024
date_stamp: "2050-01-01"
custom_sample_params: /path/to/file/with/custom_parameters.yaml
exclude_samples: [sample1, sample2]
metadata_vars: [test1, test2]
alevin:
  chemistry:
ambient:
  ambient_method: cellbender 
  cellbender_version: 'v0.3.0'
  cell_calling: barcodeRanks
  cb_max_prop_kept: 0.9 
  cb_force_expected_cells:
  cb_force_total_droplets_included:
  cb_force_low_count_threshold:
  cb_force_learning_rate:
make_sce:
  sce_bender_prob: 0.5
doublet_id:
  dbl_min_feats: 100
qc:
  qc_min_counts: 300   
  qc_min_feats: 300        
  qc_min_mito: 0             
  qc_max_mito: 0.1         
  qc_min_splice: 0          
  qc_max_splice: 1          
  qc_min_cells: 500  
integration:
  int_n_hvgs:       2000               
  int_n_dims:       50                 
  int_dbl_res:      4                   
  int_dbl_cl_prop:  0.5                 
  int_theta:        0.1                 
  int_res_ls:       [0.1, 0.2, 0.5, 1]  
  int_sel_res:      0.2
marker_genes:
  mkr_min_cl_size: 100                                  
  mkr_min_cells: 10                                      
  mkr_not_ok_re: "(lincRNA|lncRNA|pseudogene|antisense)"  
  mkr_min_cpm_mkr: 50                                     
  mkr_min_cpm_go: 1                                      
  mkr_max_zero_p: 0.5                                    
  mkr_gsea_cut: 0.1
label_celltypes:
  lbl_tissue:      "brain_cns"      
  lbl_sel_res_cl:  "RNN_snn_res.2"  
  lbl_min_pred:    0.5              
  lbl_min_cl_prop: 0.5             
  lbl_min_cl_size: 100               minimum number of cells in a cluster for that cluster to be labelled
  sce_subsets:
    oligos:       ['Oligodendrocyte', 'OPC_COP']
    microglia:    ['Micro_Mono']
    astrocytes:   ['Astrocyte', 'Ependymal_ChorPlex']
    vascular:     ['Vascular_Fibro']
    T_NK_B_cell:  ['T_NK_B_cell']
    neurons:      ['Neurons']
zoom:
  cl_grp1:
    sel_cls: [cl06]  
    zoom_res: 0.2               
    n_hvgs: 1000                
    n_dims: 20                 
    min_n_sample: 10            
    min_n_cl: 100               
    n_train: 1000
metacells: 
  subsets:  ["oligos"]
  max_cells:  [100, 400]
pseudobulk:
  celltypes:  ["all", "oligos", "microglia"]
resources:
  mb_run_alevin_fry: 8192               
  mb_save_alevin_to_h5: 8192            
  mb_run_cellbender: 32768            
  mb_get_cellbender_qc_metrics: 4096    
  mb_run_scdblfinder: 4096             
  mb_combine_scdblfinder_outputs: 8192  
  mb_make_sce_object: 16384             
  mb_run_harmony: 16384                 
  mb_run_marker_genes: 16384            
  mb_html_marker_genes: 8192            
  mb_lbl_label_celltypes: 16384         
  mb_lbl_save_subset_sces: 16384        
  mb_lbl_render_template_rmd: 4096      

```
#### Required parameters

* `proj_dir`: path to workflowr project directory (can be created with `newproj` function); has to always be an absolute path
* `fastq_dir`: path to directory with .fastq files. Fastq files need to contain sample id in the name as well as `_R1_` an `_R2_` labels for read one and read two respectievelly. If this is a relative path, {{ software_name }} will assume that's in the `proj_dir`
* `full_tag`: project label; will be added to some filenames
* `short_tag`: short project label; will be added to some filenames
* `your_name`: will appear in html outputs
* `affiliation`: will appear in html outputs
* `sample_metadata`: path to .csv file with sample metadata; If this is a relative path, {{ software_name }} will assume that's in the `proj_dir`. Spaces in column names are not allowed
* `species`: has to match one of the values in the 'genome_name' column of the `setup_parameters.csv` file which is created with the `scsetup` function; 
* `date_stamp`: will be appended to names of output files. It has to be in "YYYY-MM-DD" format
* `chemistry`: Values can be '3LT', '3v2', '3v3', '3v4', '5v1', '5v2', '5v3' and 'multiome'. 'multiome' refers exclusivelly to gene expression data generated with 10x multiome kit. {{ software_name }} currently doesn not support analysis of ATACseq multiome data. 

#### Optional parameters

* `custom_sample_params`: file with optional custom parameters for each sample (you can set custom chemistry, custom ambient and custom cellbender parameters for each sample). Example of a `custom_sample_params`file:

```yaml
sample_1:
  chemistry: 5v2
  ambient:
    knee1: 4100
    shin1: 580
    knee2: 29
    shin2: 16
sample_2: 
  chemistry: 5v2
  ambient:
    knee1: 1100
    shin1: 200
    knee2: 12
    shin2: 5
sample_3: 
  cellbender:
    total_droplets_included: 20000
```

* `exclude_samples`: list of all samples that you don't want to include in the analysis
* `metadata_vars`: list of all metadata variables in `sample_metadata` that will be used for making plots (which plots, in which reports)

##### ambient

* `ambient_method`: options are `cellbender`, `none` or `decontx` Default is `cellbender`(maybe change this)
* `cellbender_version`: options are `'v0.3.0'` and `'v0.2.0'` (maybe not necessary to have to versions available?); this parameter is only considered if `cellbender` is selected as ambient method. 
* `cb_force_expected_cells` : force `--expected-cells` parameter to be the same for all samples; only valid if `cellbender`is selected as the ambient method
* `cb_force_total_droplets_included`: force `--total-droplets-included` parameter to be the same for all samples; only valid if `cellbender`is selected as ambient method
* `cb_force_low_count_threshold`: force `--low-count-threshold` parameter to be the same for all samples; only valid if `cellbender`is selected as ambient method
* `cb_force_learning_rate`: force `--learning-rate`parameter to bet the same for all samples; only valid if `cellbender` is selected as ambient method.
* `cb_max_prop_kept`: default is 0.9; this is only relevant if ambient method is `cellbender`. This excludes samples for which cellbender calls > 90% of total_droplets_indluded as cells. 
* `cell_calling`: this param is only considered when ambient method is `none`or `decontx`. Options include `barcodeRanks` and `emptyDrops`

##### make_sce

* `sce_bender_prob`: what probability of being a cell is required to be kept as a cell? Only relevant when selected ambient method is `cellbender`.

##### doublet_id
* `dbl_min_feats`: How many features does a cell need to have to be included in scDblFinder calculations?

##### qc

* `qc_min_counts`thresholding for cells to keep: minimum number of UMIs per cell
* `qc_min_feats`thresholding for cells to keep: minimum number of features per cell
* `qc_min_mito`thresholding for cells to keep: minimum proportion of mitochondrial reads 
* `qc_max_mito`thresholding for cells to keep: maximum proportion of mitochondrial reads
* `qc_min_splice`thresholding for cells to keep: minimum proportion of spliced reads
* `qc_max_splice`thresholding for cells to keep: maximum proportion of spliced reads
* `qc_min_cells`: samples are only kept if they have at least this many cells after QC

##### integration

* `int_n_hvgs`: number of HVGs to use for PCA
* `int_n_dims`number of PCs to use for integration
* `int_dbl_res`resolution for clustering to get extra doublets
* `int_dbl_cl_prop` proportion of doublet cells in a cluster to exclude that cluster
* `int_theta` theta parameter for Harmony integration. 0 means no extra mixing of batch variable, 2 is default, which encourages batch mixing. At the moment our default is 0.1.
* `int_res_ls` list of cluster resolutions for Harmony clustering
* `int_sel_res` selected cluster resolution (for marker genes?)

##### marker_genes
* `mkr_min_cl_size` minimum number of cells in a cluster for that cluster to have marker genes calculated
* `mkr_min_cells` minimum number of cells present in a pseudobulk sample for that sample to be used to calculate marker genes
* `mkr_not_ok_re` regular expression for gene types to exclude from marker gene plotting
* `mkr_min_cpm_mkr` minimum number of counts per million in a cell type for a gene to be considered a marker gene
* `mkr_min_cpm_go` minimum number of counts per million in a cell type for a gene to be used in GO analysis
* `mkr_max_zero_p` maximum proportion of pseudobulk samples for a cell type that can have zero counts for a gene to be used in GO analysis
* `mkr_gsea_cut` FDR cutoff for GSEA analysis

##### label_celltypes

* `lbl_tissue`: which tissue do you want to label? (at present we only offer "brain_cns" :D)
* `lbl_sel_res_cl`: selected cluster resolution for cell type labelling (higher should be better)
* `lbl_min_pred`: minimum probability of a cell being a cell type for that cell to be labelled as that cell type
* `lbl_min_cl_prop`: minimum proportion of cells in a cluster that need to be labelled for that cluster to be labelled
* `lbl_min_cl_size`: minimum number of cells in a cluster for that cluster to be labelled
* `sce_subsets`: specification of sce subsets to be saved. spec should be name, followed by list of celltypes to include. valid labels for brain_cns: 'OPC_COP', 'Oligodendrocyte', 'Astrocyte', 'Ependymal', 'Choroid_plexus', 'Micro_Mono', 'Vascular_Fibro', 'T_NK_B_cell', 'Neurons'

##### zoom
* `sel_cls: [cl06]`: which clusters to include?  
* `zoom_res`: what resolution?
* `n_hvgs`: how many HVGs to use for PCA
* `n_dims`: how many PC dimensions to integrate over
* `min_n_sample`: what is the minimum number of cells needed to keep a sample?
* `min_n_cl`:  what is the minimum number of cells needed to keep a cluster?
* `n_train`: how many cells per cluster should we use for training

##### metacells
* `subset`: 
* `max_cells`:

##### pseudobulk
* `celltypes`:


##### Resources
* `retries`
* `mb_run_alevin_fry`: values are in MB for the whole job, not per thread
* `mb_save_alevin_to_h5`
* `mb_run_ambient`
* `mb_get_ambient_qc_metrics`
* `mb_run_scdblfinder`
* `mb_combine_scdblfinder_outputs`
* `mb_make_sce_object`
* `mb_run_harmony`
* `mb_run_marker_genes`
* `mb_lbl_label_celltypes`
* `mb_lbl_save_subset_sces`
* `mb_lbl_render_template_rmd`


