# Getting started


## Prerequisites

### Hardware

{{ software_name }} is intended for Linux-based high-performance computing clusters (HPCs). It will also run on workstations (i.e. powerful local machines).

If you want to use `cellbender` for ambient RNA correction, then you will need to have access to GPUs with CUDA support on your cluster. If you don't have access but still want to do ambient RNA correction, you can use `decontx`.

### Software

!!! Warning "Ignore this if on Roche shpc"

{{ software_name }} requires `snakemake` and `conda`. See the [snakemake manual](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and the [conda user guide](https://docs.anaconda.com/miniconda/) for help with the installation.

## Installation


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
    #     qos: '1d'
    #     runtime: 720
    #   run_cellbender:
    #     slurm_partition: "batch_gpu"
    #     slurm_extra: "'--gpus=1'"
    #     qos: '1d'
    #     runtime: 720
    ```

2.  Set up a conda environment that you will use for running {{software_name}}. In the {{software_name}} folder, there is a folder called _envs_, and within that you can choose between envs for `slurm`, `lsf`, and `local`.

    ```bash
    ml Miniforge3
    conda env create -n scprocess -f envs/scprocess_slurm.yaml
    # alternatives:
    # conda env create -n scprocess -f envs/scprocess_lsf.yaml
    # conda env create -n scprocess -f envs/scprocess_local.yaml
    ```

    Check that this worked:

    ```bash
    # reload the .bashrc file
    conda activate scprocess
    ```


3.  Add {{ software_name }} to your path. Open your `.bashrc` file and add the following line:

    ```bash
    export PATH=/PATH/TO/YOUR/FOLDER/scprocess:${PATH}
    ```

    Check that this worked:

    ```bash
    # reload the .bashrc file
    source ~/.bashrc

    # check that scprocess and scsetup work
    scprocess -h
    scsetup -h
    ```


### General installation instructions

1. Clone the github repository:

    ```
    git clone https://github.com/marusakod/scprocess_test.git
    ```

2.  Set up a conda environment that you will use for running {{software_name}}. In the {{software_name}} folder, there is a folder called _envs_, and within that you can choose between envs for `slurm`, `lsf`, and `local`.

    ```bash
    conda env create -n scprocess -f envs/scprocess_slurm.yaml
    # alternatives:
    # conda env create -n scprocess -f envs/scprocess_lsf.yaml
    # conda env create -n scprocess -f envs/scprocess_local.yaml
    ```

    Check that this worked:

    ```bash
    # reload the .bashrc file
    conda activate scprocess
    ```

3.  Add {{ software_name }} to your path. Open your `.bashrc` file and add the following line:

    ```bash
    export PATH=/PATH/TO/YOUR/FOLDER/scprocess:${PATH}
    ```

    Check that this worked:

    ```bash
    # reload the .bashrc file
    source ~/.bashrc

    # check that scprocess and scsetup work
    scprocess -h
    scsetup -h
    ```


## {{ software_name }} data directory setup

1. Create a directory that will store all necessary data for running {{ software_name }}:

    ```bash
    mkdir /path/to/scdata
    ```


2. Add the following line to your `.bashrc` file:

    ```bash
    export SCPROCESS_DATA_DIR=/path/to/scdata
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

!!! warning "Roche sHPC only"
    On the Roche sHPC, you should use the profile _slurm_shpc_.

{{ software_name }} is intended to be used with a cluster with a job scheduler such as `slurm` or `LSF` (although it will still work without a job scheduler). To set up a job scheduler in `snakemake`, it is common to define a configuration profile with cluster settings e.g. resource allocation. {{ software_name }} comes with two predefined configuration profiles stored in the _profiles_ directory: _profiles/slurm_default_ and _profiles/lsf_default_ for `slurm` and `LSF` respectively. 

To use {{ software_name }} with a job scheduler, you need to add a line to your  _scprocess_setup.yaml_ file:

=== slurm
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
