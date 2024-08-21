# Getting started

## Prerequisits
### System requirements
#### Hardware
scprocess runs on Linux systems that meet these minimum requirements:
* processor?
* RAM? 
* OS? Linux
* GPU with CUDA support (only required if you select CellBender as ambient method)

This is what they have for Cell Ranger:
Cell Ranger ARC pipelines run on Linux systems that meet these minimum requirements:

* 8-core Intel or AMD processor (48 cores recommended)
* 64GB RAM (320GB recommended)
* 1TB free disk space
* 64-bit CentOS/RedHat 6.0 or Ubuntu 12.04; See the 10x Genomics OS Support page for details.

#### Software

scprocess requires snakemake and conda. See the [snakemake manual](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and the [conda user guide](https://docs.anaconda.com/miniconda/) for help with the installation.

## Installation

Clone the github repository:
```
git clone https://github.com/wmacnair/scprocess.git

```
Add scprocess to your path. Go into your `.bashrc` file and add the following line:

```bash
export PATH=/PATH/TO/YOUR/FOLDER/scprocess:${PATH}
```
## Scprocess data directory setup

Create a directory that will store all data necessary for running scprocess. Add a path to that directory to your `.bashrc` file using:

```bash

export SCPROCESS_DATA_DIR=/path/to/scprocess_data_directory

```
Create a config file. Example file:

```yaml
genome:
  tenx:
    name: human_2024 
    decoys: True 
  custom:
    name: 'mouse' 
    fasta: '/path/to/fasta/genome.fa'
    gtf: '/path/to/genes/genes.gtf'
    decoys: False 
    mito_str: "^mt-" 
```

Valid values for tenx genome names are: `human_2024`, `human_2020`, `mouse_2024` and `mouse_2020`. Multiple names can be specified; duplicated names are not allowed! Describe what tenx means and what custom means (where do tenx genomes come from)

Specifying `decoys` is optional and will default to `True` for all genomes.

Add a note here, referencing something from Rob showing that decoys make sense.

Example config file for multiple genomes:

```yaml
genome:
  tenx:
    name: [human_2024, mouse_2024] 
    decoys: [True, False]
  custom:
    name: ['mouse', 'zebrafish']
    fasta: ['/path/to/mouse/fasta/genome.fa', 'path/to/zebrafish/fasta/genome.fa']
    gtf: ['/path/to/mouse/genes/genes.gtf', 'path/to/zebrafish/genes/genes.gtf']
    decoys: [True, False] 
    mito_str: ["^mt-"]
```

Fininsh setting up scprocess data directory using:

```
scsetup /path/to/setup_config.yaml

```

## Project directory setup


scprocess relies on the `workflowr` project directory structure to ensure that analyses are systematically organized and sherable. [Workflowr](https://workflowr.github.io/workflowr/) is an R package that enhances project management by integrating R markdown for literate programming, Git for version control, and automated website generation for presenting results. By incorporating workflowr into scprocess, users benefit from a streamlined workflow where each analysis step is documented, versioned and linked to it's coresponding code and environment.

```
newproj project_name -w /path/to/project/directory

```


