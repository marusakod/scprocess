# Getting started


## Prerequisites

### Hardware

{{ software_name }} runs on Linux systems that meet these minimum requirements:
 
* Operating system: Linux
* [processor? List all processors that we tested on? The only problematic part is probably alevin]
* [RAM? depends on how big the dataset is]
* [CPU?]
* [GPU with CUDA support (only required if you select CellBender as ambient method)]

### Software

!!! Warning "Ignore this if on Roche shpc"

{{ software_name }} requires `snakemake` and `conda`. See the [snakemake manual](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and the [conda user guide](https://docs.anaconda.com/miniconda/) for help with the installation.

## Installation

### General installation instructions

1. Clone the github repository:

    ```
    git clone https://github.com/marusakod/scprocess_test.git
    ```

2.  Add {{ software_name }} to your path. Open your `.bashrc` file and add the following line:

    ```bash
    export PATH=/PATH/TO/YOUR/FOLDER/scprocess:${PATH}
    ```


### Roche sHPC installation

1. Clone the Roche GitLab repository:

    ```bash
    cd ~/packages/ # or wherever you keep your packages
    git clone https://code.roche.com/macnairw/scprocess
    ```

    You should be able to see a `slurm_shpc` folder in the `profiles` folder in the top level of the {{ software_name }} directory:

    ```bash
    cat profiles/slurm_shpc/config.yaml
    # # General configurations
    # executor: slurm
    # latency-wait: 10
    # show-failed-logs: True
    # keep-going: True
    # scheduler: greedy
    # printshellcmds: True
    # jobs: 20
    # default-resources:
    #   runtime: 3h
    #   mem_mb: 4096
    # set-resources:
    #   run_harmony:
    #     qos: 1d
    #     runtime: 12h
    #   run_ambient:
    #     slurm_extra: "'-p batch_gpu --gpus=1'"
    #   run_cellbender:
    #     slurm_extra: "'-p batch_gpu --gpus=1'"
    ```

2. Add some things to your `~/.bashrc` (this code adds some extra lines to the end of your `.bashrc` file. Feel free to put them somewhere more tidy!):

    ```bash
    # add scprocess to path
    echo "export PATH=~/packages/scprocess:${PATH}" >> ~/.bashrc

    # add some sHPC-specific things
    echo "alias scprocess='export ROCS_ARCH=sandybridge; source /apps/rocs/init.sh; ml snakemake-slurm/0.15.0-foss-2020a-Python-3.11.3-snakemake-8.30.0; scprocess'" >> ~/.bashrc
    echo "alias scsetup='export ROCS_ARCH=sandybridge; source /apps/rocs/init.sh; ml snakemake-slurm/0.15.0-foss-2020a-Python-3.11.3-snakemake-8.30.0; scsetup'" >> ~/.bashrc

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

3. Create a configuration file _scprocess_setup.yaml_ in the `$SCPROCESS_DATA_DIR` directory you just created, with the contents as follows:

    ```yaml
    genome:
      tenx:
      - name: human_2024 
      - name: mouse_2024 
    ```

    This will ask the setup process to download and prepare the most recent pre-built [human](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads:~:text=Human%20reference%20(GRCh38)%20%2D%202024%2DA) and [mouse](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads:~:text=Mouse%20reference%20(GRCm39)%20%2D%202024%2DA) genomes from 10x Genomics. For more information on how to structure the _scprocess_setup.yaml_ see the [`Reference`](reference.md#scsetup) section.

    ??? tip "Save some space by removing the reference genome used for the tutorial"
        [Quick start tutorial](tutorial.md) section demonstrates how to run {{ software_name }} on an example human dataset. In order for users to be able to follow the tutorial `scsetup` will automatically download the `human_2024` reference genome and build an alevin index with [decoys](reference.md#scsetup) even if `human_2024` is not listed in the _scprocess_setup.yaml_ file. If you would like to remove this reference genome (after running the tutorial) use:
    
        ```bash
        rm -rf $SCPROCESS_DATA_DIR/reference_genomes/human_2024
        rm -rf $SCPROCESS_DATA_DIR/alevin_fry_home/human_2024

        # optionally modify the setup_parameters.csv file where all genomes available in $SCPROCESS_DATA_DIR are listed
    
        awk -F',' '$1 != "human_2024"' $SCPROCESS_DATA_DIR/setup_parameters.csv > $SCPROCESS_DATA_DIR/temp.csv && mv $SCPROCESS_DATA_DIR/temp.csv $SCPROCESS_DATA_DIR/setup_parameters.csv 
        ```
    

4. Finish setting up scprocess data directory with:

    We find it good practice to first do a "dry run" to check what will happen:
    ```bash
    scsetup -n
    # or
    # scsetup --dry-run
    ```

    If that looks ok, then run it for real:
    ```bash
    scsetup
    ```


## Cluster setup

{{ software_name }} is intended to be used with a cluster with a job scheduler such as `Slurm` or `LSF` (although it will still work without a job scheduler). To set up a job scheduler in `snakemake`, it is common to define a configuration profile with cluster settings e.g. resource allocation. {{ software_name }} comes with two predefined configuration profiles stored in the _profiles_ directory: _profiles/slurm_default_ and _profiles/lsf_default_ for `Slurm` and `LSF` respectively. 

To use {{ software_name }} with a job scheduler, you need to add a line to your  _scprocess_setup.yaml_ file:

=== Slurm
```yaml
profile: slurm_default
```
=== LSF
```yaml
profile: lsf_default
```

If you want to make a profile that is specific to your own cluster, we recommend that you make a copy one of the default profile folders, e.g. to _profiles/slurm_my_cluster_, then edit the corresponding _config.yaml_ file. Once you are happy with it, edit the _scprocess_setup.yaml_ file to point to this profile like before, e.g. 

```yaml
profile: slurm_my_cluster
```

`scsetup` and {{ software_name }} will then run in cluster mode with the specifications in this profile.
