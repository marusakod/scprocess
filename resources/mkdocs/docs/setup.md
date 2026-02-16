# Getting started


## Prerequisites

### Hardware

{{sc}} is designed to run on Linux-based high-performance computing (HPC) clusters but can also be used on powerful standalone workstations.

For users intending to perform ambient RNA correction with `cellbender`, access to GPU with CUDA support is required. If you do not have access to suitable GPU but still wish to perform ambient RNA correction, you can use the alternative tool, `decontx`.


### Software

{{sc}} requires Conda. If you do not already have Conda installed, refer to the [Conda user guide](https://docs.anaconda.com/miniconda/) for detailed installation instructions.

If you plan to use `cellbender` for ambient RNA correction, you will also need Apptainer. For guidance on installing Apptainer, see the [installation instructions](https://apptainer.org/docs/admin/main/installation.html).


## Installation

1. Clone the repository:

    ```bash
     git clone git@github.com:marusakod/scprocess.git
    ```
  
2.  Create a Conda environment:
    
    Navigate into the {{sc}} directory and create a Conda environment. Choose the appropriate command based on your operating environment (local machine, SLURM, or LSF cluster).
    
    === "local"
        ```bash
        conda env create -n scprocess -f envs/scprocess_local.yaml
        ```

    === "SLURM"
        ```bash
        conda env create -n scprocess -f envs/scprocess_slurm.yaml
        ```

    === "LSF"
        ```bash
        conda env create -n scprocess -f envs/scprocess_lsf.yaml
        ```

    If you are using LSF or SLURM, remember to also review the [Cluster setup](#cluster-setup) section below.

3.  Add {{sc}} to your PATH. 
    
    Open your `~/.bashrc` file and add the following line:

    ```bash
    export PATH=/PATH/TO/YOUR/FOLDER/scprocess:${PATH}
    ```
    
    Verify the installation by:

    * Reloading your `~/.bashrc` file: 
      ```bash
      source ~/.bashrc
      ```
    
    * Activating the {{sc}} Conda environment and checking for help messages:
      ```bash
      conda activate scprocess
      scprocess -h
      ```


## {{sc}} data directory setup

{{sc}} requires a dedicated directory to store all necessary data, such as reference genomes.

1. Create the data directory:

    ```bash
    mkdir /path/to/scdata
    ```

2. Add the following line to your `.bashrc` file:

    ```bash
    export SCPROCESS_DATA_DIR=/path/to/scdata
    ```

3. Create a configuration file _scprocess_setup.yaml_ in the `$SCPROCESS_DATA_DIR` directory you just created, with the contents as follows:

    ```yaml
    user:
      your_name:    Testy McUser
      affiliation:  Unemployed
    ref_txomes:
      tenx:
      - name: human_2024 
      - name: mouse_2024 
    ```

    This will ask the setup process to download and prepare the most recent pre-built [human](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads:~:text=Human%20reference%20(GRCh38)%20%2D%202024%2DA) and [mouse](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads:~:text=Mouse%20reference%20(GRCm39)%20%2D%202024%2DA) reference transcriptomes from 10x Genomics. For more information on how to structure the _scprocess_setup.yaml_ see the [`Reference`](reference.md#scprocess-setup) section.

    ??? tip "Save some space by removing the reference transcriptome used for the tutorial"
        [Quick start tutorial](tutorial.md) section demonstrates how to run {{sc}} on an example mouse datasets. In order for users to be able to follow the tutorial `scprocess setup` will automatically download the `mouse_2024` reference transcriptome. If you would like to remove it (after running the tutorial) use:
    
        ```bash
        rm -rf $SCPROCESS_DATA_DIR/reference_transcriptomes/mouse_2024
        rm -rf $SCPROCESS_DATA_DIR/alevin_fry_home/mouse_2024

        # optionally modify the setup_parameters.csv file where all genomes available in $SCPROCESS_DATA_DIR are listed
    
        awk -F',' '$1 != "mouse_2024"' $SCPROCESS_DATA_DIR/setup_parameters.csv > $SCPROCESS_DATA_DIR/temp.csv && mv $SCPROCESS_DATA_DIR/temp.csv $SCPROCESS_DATA_DIR/index_parameters.csv 
        ```
    

4. Finish setting up the {{sc}} data directory:

    To download all required data and index reference transciptomed use the {{scsetup}} command. The first time you run {{scsetup}} you need to specify a `-c`/`--rangerurl` flag and provide a valid download link for Cellranger (v9.0.0 or higher) available on the [10x Genomics CellRanger download & installation page](https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions) : 
    
    ```bash
    scprocess setup -c "https://cf.10xgenomics.com/releases/cell-exp/cellranger-10.0.0.tar.gz..." 
    ```
    Once the inital setup is complete, you do not need to provide the CellRanger link again i.e if you modify the _scprocess_setup.yaml_ file to add additional reference genomes, simply run {{scsetup}}.

## Cluster setup

{{sc}} is intended to be used on a cluster with a job scheduler such as `SLURM` or `LSF` (although it will also work without a job scheduler). To set up a job scheduler in `snakemake`, it is common to define a configuration profile with cluster settings e.g. resource allocation. {{sc}} comes with two predefined configuration profiles stored in the _profiles_ directory: _profiles/slurm_default_ and _profiles/lsf_default_ for `SLURM` and `LSF` respectively. 

To use {{sc}} with a job scheduler, you need to add a line to your  _scprocess_setup.yaml_ file:

=== "SLURM"
```yaml
user:
  profile: slurm_default
```
=== "LSF"
```yaml
user:
  profile: lsf_default
```

If you want to make a profile that is specific to your cluster, we recommend that you make a copy one of the default profile folders, e.g. to _profiles/slurm_my_cluster_, then edit the corresponding _config.yaml_ file. Once you are happy with it, edit the _scprocess_setup.yaml_ file to point to this profile like before, e.g. 

```yaml
user:
  profile: slurm_my_cluster
```
{{scsetup}} and {{scrun}} will then run in cluster mode with the specifications in this profile.
