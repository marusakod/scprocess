# Running scprocess

To run scprocess do:
```
scprocess /path/to/config.yaml
```

If you want to run a dry run you can add a `-n` or `--dry-run` flag to this command.
In case you need to use other snakemake options that are not included in scprocess, you can use the `-E` or `--extraargs` flag. It accepts any text and converts it into a string. The argument of `-E` has to be between quotes. For example, if you need to run a dry-run (use another example cause there is already -n for dryrun): 

```
scprocess $CONFIGFILE -E " -np "
```
Add here that it's useful to add a max memory and max threads argument here (If you're not using a cluster). `--cores` will overwrite all threads specified in rules (what if we want to everwrite just what is > than max_cores?)

### config file

#### Required parameters

This is an example config file for scprocess with minimal required parameters:

```yaml

proj_dir: /path/to/proj/directory
fastq_dir: /path/to/directory/with/fastq/files
full_tag: test_project_204
short_tag: test
your_name: John Doe
affiliation: where you work
sample_metadata: /path/to/metadata
alevin:
  species: human_2024
```

* `proj_dir`: path to workflowr project directory (can be created with `newproj` function); has to always be an absolute path
* `fastq_dir`: path to directory with .fastq files. Fastq files need to contain sample id in the name as well as '_R1_' an '_R2_' labels for read one and read two respectievelly. If this is a relative path, scprocess will assume that's in the `proj_dir`
* `full_tag`: Project label; will be added to some filenames
* `short_tag`: short project label; will be added to some filenames
* `your_name`: will appear in html outputs
* `affiliation`: will appear in html outputs
* `sample_metadata`: path to .csv file with sample metadata; has to contain column named 'sample_id' with values matching sample ids in fastq files; If this is a relative path, scprocess will assume that's in the `proj_dir`
* `species`: has to match one of the values in the 'genome_name' column of the `setup_parameters.csv` file which is created with the `scsetup` function; (it might make sense to move species outside of the alevin section, next to other reqired params?)

#### Optional parameters

```yaml

date_stamp: "2050-01-01"
exclude_samples: [sample1, sample2]
metadata_vars: [test1, test2]
```

* `date_stamp` you can add a specific date to the config which will be appended to output filenames. If you don't specify this in the config, scprocess will use the currect date (which is annoying if it runs for more than 1 day cause than it will start from scratch.)
* `exclude_samples`: list of all samples that you don't want to include in the analysis
* `metadata_vars`: list of all metadata variables in `sample_metadata` that will be used for making plots (which plots? in which .html report?)


##### Simpleaf

```yaml

alevin:
  species: human 
  chemistry: /path/to/sample/chemistry_file.csv
```
scprocess will detect the chemistry automatically. Alternatively you can add a `chemistry` parameter to the config file which should be path to csv file with columns 'sample_id' and 'version'; valid versions are '3LT', '3v2', '5v1', '5v2', '3v3', 'multiome'. Add a note that 'multiome' here is referring exclusivelly to gene expression data generated with 10x multiome kit, and that scprocess currently doesn not support analysis of ATACseq multiome data. 


##### Ambient RNA removal

```yaml

ambient:
  ambient_method: cellbender 
  cellbender_version: 'v0.3.0'
  custom_params: /path/to/file/with/custom_parameters.csv
  cb_max_prop_kept: 0.9 
  cell_calling: barcodeRanks 

```

* `ambient_method`: options are `cellbender`, `none` or `decontx` Default is `cellbender`(maybe change this)
* `cellbender_version`: options are `'v0.3.0'` and `'v0.2.0'`; (check if those are really version names; Why do we have 2 versions available); this parameter is only considered if `cellbender` is selected as ambient method. 
* `custom_params`: describe what the format of this file should be like for each individual method (add examples and show on one example knee plot how the thresholds should be set)

if ambient method is `cellbender`, the `custom_params` file should include columns 'sample_id', 'total_droplets_included', 'expected_cells', 'low_count_threshold', 'learning_rate', 'empty_start' and 'empty_end'. (empty start and empty end are actually needed for the pb_empties rule not here)

if another ambient method is selected than the format of the `custom_params` file depends on the method selected for cell calling. If `cell_calling` is set to `barcodeRanks`, column names need to include 'expected_cells', 'empty_start' and 'empty_end'. If `cell_calling` is set to `emptyDrops`, columna names need to include 'retain', 'empty_start' and 'empty_end'. 

* `cb_max_prop_kept`: default is 0.9; this is only relevant if ambient method is `cellbender`. This excludes samples for which cellbender calls > 90% of total_droplets_indluded as cells. 
* `cell_calling`: this param is only considered when ambient method is `none`or `decontx`. Options include `barcodeRanks` and `emptyDrops`


##### Make sce object

```yaml

make_sce:
  sce_bender_prob: 0.5

```
you can set the `sce_bender_prob` parameter which is only relevant when selected ambient method is cellbender. Meaning: what probability of being a cell is required to be kept as a cell?

##### Doublet removal

```yaml 
doublet_id:
  dbl_min_feats = 100 
```

* `dbl_min_feats`: How many features does a cell need to have to be included in scDblFinder calculations?

##### QC
```yaml
qc:
  qc_min_counts: 300   
  qc_min_feats: 300        
  qc_min_mito: 0             
  qc_max_mito: 0.1         
  qc_min_splice: 0          
  qc_max_splice: 1          
  qc_min_cells: 500        

```

* `qc_min_counts`thresholding for cells to keep: minimum number of UMIs per cell
* `qc_min_feats`thresholding for cells to keep: minimum number of features per cell
* `qc_min_mito`thresholding for cells to keep: minimum proportion of mitochondrial reads 
* `qc_max_mito`thresholding for cells to keep: maximum proportion of mitochondrial reads
* `qc_min_splice`thresholding for cells to keep: minimum proportion of spliced reads
* `qc_max_splice`thresholding for cells to keep: maximum proportion of spliced reads
* `qc_min_cells`: samples are only kept if they have at least this many cells after QC

##### Integration

```yaml

integration:
  int_n_hvgs:       2000               
  int_n_dims:       50                 
  int_dbl_res:      4                   
  int_dbl_cl_prop:  0.5                 
  int_theta:        0.1                 
  int_res_ls:       [0.1, 0.2, 0.5, 1]  
  int_sel_res:      0.2                 

```

* `int_n_hvgs`: number of HVGs to use for PCA
* `int_n_dims`number of PCs to use for integration
* `int_dbl_res`resolution for clustering to get extra doublets
* `int_dbl_cl_prop`proportion of doublet cells in doublet cluster to exclude that cluster
* `int_theta` theta parameter for Harmony integration. 0 means no extra mixing of batch variable, 2 is default, which encourages batch mixing. At the moment our default is 0.1.
* `int_res_ls` list of cluster resolutions for Harmony clustering
* `int_sel_res` selected cluster resolution (for marker genes?)


##### Marker genes

```yaml

marker_genes:
  mkr_min_cl_size: 100                                  
  mkr_min_cells: 10                                      
  mkr_not_ok_re: "(lincRNA|lncRNA|pseudogene|antisense)"  
  mkr_min_cpm_mkr: 50                                     
  mkr_min_cpm_go: 1                                      
  mkr_max_zero_p: 0.5                                    
  mkr_gsea_cut: 0.1                                       

```

* `mkr_min_cl_size` minimum number of cells in a cluster for that cluster to have marker genes calculated
* `mkr_min_cells` minimum number of cells present in a pseudobulk sample for that sample to be used to calculate marker genes
* `mkr_not_ok_re` regular expression for gene types to exclude from marker gene plotting
* `mkr_min_cpm_mkr` minimum number of counts per million in a cell type for a gene to be considered a marker gene
* `mkr_min_cpm_go` minimum number of counts per million in a cell type for a gene to be used in GO analysis
* `mkr_max_zero_p` maximum proportion of pseudobulk samples for a cell type that can have zero counts for a gene to be used in GO analysis
* `mkr_gsea_cut` FDR cutoff for GSEA analysis


##### Label celltypes

```yaml

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

```

* `lbl_tissue`: which tissue do you want to label? (at present we only offer "brain_cns" :D)
* `lbl_sel_res_cl`: selected cluster resolution for cell type labelling (higher should be better)
* `lbl_min_pred`: minimum probability of a cell being a cell type for that cell to be labelled as that cell type
* `lbl_min_cl_prop`: minimum proportion of cells in a cluster that need to be labelled for that cluster to be labelled
* `lbl_min_cl_size`: minimum number of cells in a cluster for that cluster to be labelled
* `sce_subsets`: specification of sce subsets to be saved. spec should be name, followed by list of celltypes to include. valid labels for brain_cns: 'OPC_COP', 'Oligodendrocyte', 'Astrocyte', 'Ependymal', 'Choroid_plexus', 'Micro_Mono', 'Vascular_Fibro', 'T_NK_B_cell', 'Neurons'

##### Zoom

```yaml
zoom:
  cl_grp1:
    sel_cls: [cl06]  
    zoom_res: 0.2               
    n_hvgs: 1000                
    n_dims: 20                 
    min_n_sample: 10            
    min_n_cl: 100               
    n_train: 1000               

```

* `sel_cls: [cl06]`: which clusters to include?  
* `zoom_res`: what resolution?
* `n_hvgs`: how many HVGs to use for PCA
* `n_dims`: how many PC dimensions to integrate over
* `min_n_sample`: what is the minimum number of cells needed to keep a sample?
* `min_n_cl`:  what is the minimum number of cells needed to keep a cluster?
* `n_train`: how many cells per cluster should we use for training

##### Resources

```yaml

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

Values are in MB, for the whole job (not per thread). These are inital values. If the job fails, this value will be multiplied by 2,3,4 then 5. Then scprocess will give up (does this make sense if not run in cluster mode?)

