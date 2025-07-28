# Usage

## Basic usage

Assuming the required hardware is available, all software is installed and you have successfully completed the [setup of scprocess data directory](setup.md#scprocess-data-directory-setup), you can run {{sc}} on your dataset by following the steps outlined below:

1. [Prepare project directory](usage.md#1-prepare-project-directory)
2. [Prepare input files](usage.md#2-prepare-input-files)
3. [Prepare configuration file (config.yaml)](usage.md#3-prepare-configuration-file-configyaml)
4. [Run the analysis](usage.md#4-run-the-analysis)

!!! warning "{{sc}} expects multiple samples"

### 1. Prepare project directory

{{sc}} relies on the [`workflowr`](https://workflowr.github.io/workflowr/) project directory template. You can create a new `workflowr` project using {{scnew}}, as follows:

```bash
# create a new project in the current directory, with directories for fastq and metadata files, and a default config file
scprocess newproj my_project -c -s
```

### 2. Prepare input files

{{sc}} requires 2 types of input files:

* **FASTQ files** generated using 10x Genomics technology: names of FASTQ files have to contain a `[SAMPLE_ID]` in the name as well as `_R1_` an `_R2_` labels for read one (forward read) and read two (reverse read) respectively. For example:

    `[SAMPLE_ID]*_R1_*.fastq.gz` and `[SAMPLE_ID]*_R2_*.fastq.gz`, where `*` can be replaced with any character.

* **metadata**: a CSV file with sample information. The only required column in the metadata file is `sample_id`, where the values should match the `[SAMPLE_ID]` labels included in FASTQ file names.

To see examples of input files go to [Quick start tutorial](tutorial.md#creating-a-new-project-directory-and-preparing-input-data). 

### 3. Prepare configuration file

{{sc}} requires a configuration YAML file where you can specify analysis parameters. This is an example of a configuration file with all [required parameters](reference.md#required-parameters).

```yaml
proj_dir: /path/to/proj/directory
fastq_dir: /path/to/directory/with/fastq/files
full_tag: test_project
short_tag: test
your_name: John Doe
affiliation: where you work
sample_metadata: /path/to/metadata.csv
species: human_2024
date_stamp: "2025-01-01"
alevin:
 chemistry: 3v3
```

### 4. Run the analysis

To run {{sc}} do:

```bash
scprocess run /path/to/config-my_project.yaml
```
If you want to run a dry run you can add a `-n` or `--dry-run` flag to this command. In case you need to use other [snakemake options](https://snakemake.readthedocs.io/en/stable/executing/cli.html) that are not included in {{scrun}} by default, you can use the `-E` or `--extraargs` flag. For example, if you would like to set a global maximum for the number of threads available to any rule you can use: 

```bash
scprocess run /path/to/config-my_project.yaml -E " --max-threads 8 "
```

By default {{scrun}} will run rule `all` which includes all [core steps](introduction.md#rule-all-steps). The [optional steps](introduction.md#optional-steps) (rule `label_celltypes`) can run only after  rule `all` is completed and have to be specifically requested.

Additionally, you can run individual rules that generate HTML outputs (`mapping`, `ambient`, `qc`, `integration`, `marker_genes`). This is useful if you want to inspect the html outputs for the intermediate steps first and then continue with the analysis. To run each rule separately you have to specify the rule using the `-r` or `--rule` flag e.g.

```bash
scprocess run /path/to/config.yaml -r qc
```

## Analysis of multiplexed samples

{{sc}} supports two approaches for handling multiplexed samples:

* **Hashtag oligo (HTO)-based demultiplexing**: {{sc}} uses HTO-derived cDNA libraries to generate a count matrix which can be used for sample demultiplexing.

* **Outputs of external demultiplexing algorithms**: If the data has already been demultiplexed using an external method (e.g. genetic-based tools like [demuxlet](https://www.nature.com/articles/nbt.4042)), users can provide a cell-sample assignment file to process the data further using {{sc}}

### Input files 

Processing multiplexed samples requires a different format for the sample metadata CSV file. In addition to the standard `sample_id` column, the following columns must be included:

* `pool_id`: specifies the pool to which each sample belongs. Instead of values in the `sample_id` column, FASTQ filenames must match values in the `pool_id` column.
* `hto_id`: only required for HTO-based demultiplexing. Specifies the HTO label used to tag each sample before pooling.


!!! Warning "The dataset must consist entirely of either multiplexed or non-multiplexed samples. Mixed datasets are not supported."


![multiplexing](assets/images/scprocess_multiplexing_white_bg.png#only-light)
![multiplexing](assets/images/scprocess_multiplexing_black_bg.png#only-dark)

---

<div class="img-caption">Schematic representation of sample multiplexing for single-cell sequencing. Individual samples (with corresponding names in the <code>sample_id</code> column) are labelled with antibodies carrying different HTOs (with corresponding labels in the <code>hto_id</code> column). These labeled samples are then combined into pools (with corresponding names in the <code>pool_id</code> column). HTO labels can be shared across different pools. </div>

## Best practices

The default parameters in the configuration file are suitable for running {{sc}} on the [example dataset](tutorial.md). Here are some of the most important things to consider when you are running {{sc}} on your own dataset:

### Setting parameters for ambient contamination removal

#### Ambient method

By default {{sc}} will use `decontx` for ambient RNA removal, which doesn't require GPU. If a GPU is available, we recommend using `cellbender` for ambient RNA decontamination as it was found to perform better than other related algorithms in [the latest benchmark](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02978-x).

#### Knee parameters

![empties_cells](assets/images/knee_plot_with_cells_and_empties.png)
Both algorithms for ambient RNA decontamination available in {{sc}} estimate background noise from empty droplets. Therefore, correctly identifying the subset of barcodes corresponding to empty droplets is critical. In the barcode-rank "knee plot", where barcodes are ranked in descending order based on their library size, two distinct plateaus are typically observed: the first plateau represents droplets containing cells with high RNA content, while the second corresponds to empty droplets containing ambient RNA.

{{sc}} identifies the cell-containing and empty droplet populations by detecting key transition points on the barcode-rank curve — namely, the inflection and knee points. These points allow {{sc}} to infer the optimal parameters for both `decontx` and `cellbender`. Additionally, {{sc}} uses these estimates to identify genes enriched in empty droplets.

We recommend verifying the accuracy of these parameters by inspecting knee plots after running `mapping`. The two main parameters inferred by {{sc}} based on transition points in the barcode-rank curve are `expected_cells` and the `empty_plateau_middle` (which corresponds to the `--total-droplets-included` parameter in `cellbender`). The `empty_plateau_middle` should extend a few thousand barcodes into the second plateau.

![all_knee_examples](assets/images/all_knee_examples.png)

---

<div class="img-caption">Three examples of knee plots with four transition points (<code>knee1</code>, <code>shin1</code>, <code>knee2</code>, <code>shin2</code>) represented by purple horizontal lines and two inferred parameters (<code>expected_cells</code> and <code>empty_plateau_middle</code>) represented by blue vertical lines. In example (A) a knee plot with two clearly distinguished plateaus and properly predicted parameters is shown. In example (B) the same knee plot is shown, however the predicted parameters are wrong. Example (C) represents a sample with very high ambient RNA contamination making it impossible to distinguish the cells and empty droplets populations by eye. </div>

To identify problematic samples, {{sc}} computes two diagnostic ratios:

* `expected_cells`/`empty_plateau_middle` ratio: this helps assess whether the estimated number of cells is reasonable compared to the `empty_plateau_middle`. In example B this ratio would be increased [but so would be the slope ratio so maybe not the best example]
* `slope_ratio`: This is the ratio of the slope of the barcode-rank curve in the empty droplet region compared to the slope at the first inflection point. Samples with a high slope ratio, as seen in example C, are likely problematic because the empty droplet plateau is not clearly distinguishable. In such cases, ambient RNA contamination algorithms like `decontx` and `cellbender` may struggle to accurately estimate background noise, and we recommend considering removing these samples from further analysis.

If {{sc}} fails to estimate the knee plot parameters but the barcode-rank curve appears normal, we suggest manually adjusting the `knee1`, `knee2`, `shin1`, and `shin2` parameters in the `custom_sample_params` file. A convenient way to fine-tune these parameters is by using the [`plotknee`](reference.md#plotknee) function in {{sc}}. This allows for easy visualization and adjustment of knee points.
    

### Setting QC parameters

Depending on whether you are using {{sc}} on single cell or single nuclei data you may want to consider adjusting the default threshold for maximum allowed spliced proportion (the `qc_max_splice` parameter) in the configuration file. [Anything else?]
    
### Managing resources





