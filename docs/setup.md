# Getting started

## Prerequisites

#### Hardware

{{ software_name }} runs on Linux systems that meet these minimum requirements:
 
* Operating system: Linux
* [processor? List all processors that we tested on? The only problematic part is probably alevin]
* [RAM? depends on how big the dataset is]
* [CPU?]
* [GPU with CUDA support (only required if you select CellBender as ambient method)]

#### Software

{{ software_name }} requires `snakemake` and `conda`. See the [snakemake manual](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and the [conda user guide](https://docs.anaconda.com/miniconda/) for help with the installation.

## Installation

1. Clone the Roche gitlab repository:

    ```bash
    cd ~/packages/ # or wherever you keep your packages
    git clone https://code.roche.com/macnairw/scprocess
    ```

    You should be able to see a file _lsf.yaml_ in the top level of the {{ software_name }} directory:

    ```bash
    cat lsf.yaml
    # app_profile:
    #     - none
    # __default__:
    #   - '-q short'
    # run_ambient:
    #   - "-q short"
    #     # - "-q long"
    #     # - "-gpu 'num=1:j_exclusive=yes'"
    # run_harmony:
    #   - "-q long"
    ```

2. Add some things to your `~/.bashrc`:

    ```bash
    # add scprocess to path
    echo "export PATH=~/packages/scprocess:${PATH}" >> ~/.bashrc
    # add some sHPC-specific things
    echo "alias scprocess='export ROCS_ARCH=sandybridge; source /apps/rocs/init.sh; ml snakemake-lsf/1.0.7-foss-2020a-Python-3.11.3-snakemake-8.23.0; scprocess'" >> ~/.bashrc
    echo "alias scsetup='export ROCS_ARCH=sandybridge; source /apps/rocs/init.sh; ml snakemake-lsf/1.0.7-foss-2020a-Python-3.11.3-snakemake-8.23.0; scsetup'" >> ~/.bashrc
    # this code adds some extra lines to the end of your .bashrc file. feel free to put them somewhere more tidy!
    ```

    Check that this worked:

    ```bash
    # reload the .bashrc file
    source ~/.bashrc

    # check that scprocess works
    scprocess -h
    ```

## {{ software_name }} data directory setup

1. Create a directory that will store all necessary data for running {{ software_name }}

2. Add the following line to your `.bashrc` file:

    ```bash
    export SCPROCESS_DATA_DIR=/path/to/scprocess_data_directory
    ```

3. Create a configuration file `.scprocess_setup.yaml` in the SCPROCESS_DATA_DIR directory you just created, with the contents as follows:

    ```yaml
    genome:
      tenx:
      - name: human_2024 
      - name: mouse_2024 
    ```

    This will ask the setup process to download and prepare the most recent pre-built [human](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads:~:text=Human%20reference%20(GRCh38)%20%2D%202024%2DA) and [mouse](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads:~:text=Mouse%20reference%20(GRCm39)%20%2D%202024%2DA) genomes from 10x Genomics. For more information on how to structure the `.scprocess_setup.yaml` see the [`Reference`](reference.md#setup-config-file) section.

    ??? tip "Save some space by removing the reference genome used for the tutorial"
        [Quick start tutorial](tutorial.md) section demonstrates how to run {{ software_name }} on an example human dataset. In order for users to be able to follow the tutorial `scsetup` will automatically download the `human_2024` reference genome and build an alevin index with [decoys](reference.md#setup-config-file) even if `human_2024` is not listed in the `.scprocess_setup.yaml` file. If you would like to remove this reference genome (after running the tutorial) use:
    
        ```bash
        rm -rf $SCPROCESS_DATA_DIR/reference_genomes/human_2024
        rm -rf $SCPROCESS_DATA_DIR/alevin_fry_home/human_2024

        # optionally modify the setup_parameters.csv file where all genomes available in $SCPROCESS_DATA_DIR are listed
    
        awk -F',' '$1 != "human_2024"' $SCPROCESS_DATA_DIR/setup_parameters.csv > $SCPROCESS_DATA_DIR/temp.csv && mv $SCPROCESS_DATA_DIR/temp.csv $SCPROCESS_DATA_DIR/setup_parameters.csv 
        ```
    

4. Finish setting up scprocess data directory with:

    ```bash
    scsetup
    ```

## Cluster setup

When running {{ software_name }} on a cluster with a job scheduler like SLURM or LSF, it is common to define a configuration profile with cluster settings e.g. resource allocation. {{ software_name }} comes with two predefined configuration profiles stored in the `profiles` directory: `profiles/slurm_default` and `profiles/lsf_default` for SLURM and LSF respectively. You can add additional profiles or edit one of the profiles that already exists in {{ software_name }}. To run `scsetup` and {{ software_name }} in cluster mode add the name of the configuration profile to `.scprocess_setup.yaml` file e.g:

```yaml
profile: slurm_default
```

Note that default configuration profiles define resource requirements for default {{ software_name }} parameters. If GPU is available and you would like to select `cellbender` for [ambient RNA removal](introduction.md#ambient-rna-removal-optional), add the highlighted section to the configuration profile:

=== "profiles/slurm_default/congif.yaml"

    ```yaml hl_lines="14 15"
    executor: slurm
    latency-wait: 10
    show-failed-logs: True
    keep-going: True
    scheduler: greedy   
    printshellcmds: True
    jobs: 20
    default-resources:
      runtime: 3h
      mem_mb: 4096
    set-resources:
      run_harmony:
        runtime: 12h
      run_ambient:
        slurm_extra: "'--gpus=1'"
    ```

=== "profiles/lsf_default/congif.yaml"

    ```yaml hl_lines="14 15"
    executor: lsf
    latency-wait: 10
    show-failed-logs: True
    keep-going: True
    scheduler: greedy   
    printshellcmds: True
    jobs: 20
    default-resources:
      runtime: 3h  # change this with queue
      mem_mb: 4096 #change this
    set-resources:
      run_harmony:
        runtime: 12h # change this with queue
      run_ambient:
        slurm_extra: "'--gpus=1'" # change this
    ```




