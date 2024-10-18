# Getting started

## Prerequisits

#### Hardware

scprocess runs on Linux systems that meet these minimum requirements:

* [processor?]
* [RAM?]
* [OS? Linux]
* [GPU with CUDA support (only required if you select CellBender as ambient method)]

#### Software

scprocess requires snakemake and conda. See the [snakemake manual](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and the [conda user guide](https://docs.anaconda.com/miniconda/) for help with the installation.

1. Check that `snakemake` is working properly by entering the command `snakemake --version`. You should see something like `8.7.0` as output.

## Installation

1. Clone the github repository:

  ```
  git clone https://github.com/wmacnair/scprocess.git
  ```

2.  Add scprocess to your path. Open your `.bashrc` file and add the following line:

  ```bash
  export PATH=/PATH/TO/YOUR/FOLDER/scprocess:${PATH}
  ```

## Scprocess data directory setup

1. Create a directory that will store all data necessary for running scprocess

2. Add the following line to your `.bashrc` file:

```bash
export SCPROCESS_DATA_DIR=/path/to/scprocess_data_directory
```

3. Create a configuration file `.scprocess_setup.yaml` in the scprocess data directory you just created, with the contents as follows:

```yaml
genome:
  tenx:
  - name: human_2024 
  - name: mouse_2024 
```

This will ask the setup process to download and prepare the most recent pre-built human and mouse genomes from 10x Genomics.

You can select any of the pre-built human or mouse reference genomes from 10x Genomics by adding their names to the `tenx` section of the configuration file. Valid names for 10x Genomics reference genomes include: `human_2024`, `human_2020`, `mouse_2024`, and `mouse_2020`.

[make sure next part has moved to #reference]

To use custom reference genomes, you need to provide the following parameters in the `custom` section of the configuration file: 

* Paths to FASTA and GTF files.
* The `mito_str` parameter, a regular expression pattern used to identify mitochondrial genes in your genome annotation.

The `decoys` parameter is optional. If not specified, it defaults to `True` for all genomes.

!!! info "More about decoys"
  scprocess utilizes simpleaf, a lightweight mapping approach that, by default, maps sequenced fragments exclusively to the transcriptome. However, this can lead to incorrect mapping of reads that arise from unannotated genomic loci to the transcriptome. To mitigate this issue, the `decoys` parameter in `scsetup` it to `True`. This option allows simpleaf to identify genomic regions with sequences similar to those in transcribed regions (decoys), thereby reducing the likelihood of false mappings. We strongly recommend keeping the decoy setting enabled. For further details, refer to [Srivastava et al., 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8).

4. Finish setting up scprocess data directory with:

```bash
scsetup
# you can do a test run to see what the setup process would do by adding -n, or --dry-run
```

## Project directory setup

scprocess relies on the ['Workflowr'](https://workflowr.github.io/workflowr/) project directory template. You can create a new `workflowr` project using:

```
newproj project_name -w /path/to/project/directory
```
