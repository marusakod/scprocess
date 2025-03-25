# Getting started

## Installation

1. Clone the Roche GitLab repository:

    ```bash
    cd ~/packages/ # or wherever you keep your packages
    git clone https://code.roche.com/macnairw/scprocess
    ```
    Swith to `main-shpc` branch:
    
    ```bash
    git checkout main-shpc
    ```

    You should be able to see a `lsf.yaml` file in the top level of the {{ software_name }} directory:

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

2. Add some things to your `~/.bashrc` (this code adds some extra lines to the end of your `.bashrc` file. Feel free to put them somewhere more tidy!):

    ```bash
    # add scprocess to path
    echo "export PATH=~/packages/scprocess:${PATH}" >> ~/.bashrc

    # add some sHPC-specific things
    echo "alias scprocess='export ROCS_ARCH=sandybridge; source /apps/rocs/init.sh; ml snakemake-lsf/1.0.7-foss-2020a-Python-3.11.3-snakemake-8.23.0; scprocess'" >> ~/.bashrc
    echo "alias scsetup='export ROCS_ARCH=sandybridge; source /apps/rocs/init.sh; ml snakemake-lsf/1.0.7-foss-2020a-Python-3.11.3-snakemake-8.23.0; scsetup'" >> ~/.bashrc

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

## Cellbender setup

Note that default configuration profile (`lsf.yaml` file) defines resource requirements for default {{ software_name }} parameters. If you would like to select `cellbender` for [ambient RNA removal](introduction.md#ambient-rna-removal-optional), make sure your `lsf.yaml` contains this highlighted section:

```yaml hl_lines="5 6 7"
app_profile:
    - none
__default__:
  - '-q short'
run_ambient:
  - "-q long"
  - "-gpu 'num=1:j_exclusive=yes'"
run_harmony:
  - "-q long"
```



