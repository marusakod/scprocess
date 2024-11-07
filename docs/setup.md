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

1. Clone the github repository:

    ```
    git clone https://github.com/marusakod/scprocess_test.git
    ```

2.  Add {{ software_name }} to your path. Open your `.bashrc` file and add the following line:

    ```bash
    export PATH=/PATH/TO/YOUR/FOLDER/scprocess:${PATH}
    ```

## {{ software_name }} data directory setup

1. Create a directory that will store all necessary data for running {{ software_name }}

2. Add the following line to your `.bashrc` file:

    ```bash
    export SCPROCESS_DATA_DIR=/path/to/scprocess_data_directory
    ```

3. Create a configuration file `.scprocess_setup.yaml` in the {{ software_name }} data directory you just created, with the contents as follows:

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

When running {{ software_name }} on a cluster with a job scheduler like SLURM or LSF, it's common to define a configuration profile. A profile defines default cluster settings, such as resource allocation (CPUs, memory, runtime) and job submission commands. {{ software_name }} comes with two predefined configuration profiles stored in the `profile` directory: `profile/slurm_default` and `profile/lsf_default` for SLURM and LSF respectively. You can add additional profiles or edit one of the profiles that already exists in {{ software_name }}. [The structure of the profile will depend on the version of `snakemake` and cluster settings.]

`scsetup` and `scprocess` commands will run in cluster mode if you add the `-E " --workflow-profile={profile_name} "` option. Example for SLURM:

```bash
scsetup -E " --workflow-profile=slurm_default "
```



=== "profile/slurm_default/congif.yaml"

    ```yaml
    set-resources:
      big_job:
        cpus_per_task: 3
        mem_mb: 2000
        runtime: "1h"
      small_job:
        cpus_per_task: 1
        mem_mb: 1000
        runtime: 10
    ```

=== "profile/lsf_default/congif.yaml"

    ```yaml
    snakefile: Snakefile
    cores: 1
    latency-wait: 10
    reason: True
    show-failed-logs: True
    keep-going: True
    printshellcmds: True
    rerun-incomplete: True
    restart-times: 3
    jobname: "{rule}.{jobid}"           
    max-jobs-per-second: 1                
    max-status-checks-per-second: 10      
    jobs: 100                              
    cluster: "bsub --output=\"{proj_dir}/.log/{rule}/sample={wildcards.sample}/%J.out\" --error=\"{proj_dir}/.log/{rule}/sample={wildcards.sample}/%J.err\""

    ```




