
# Setting up `scprocess` on sHPC

## First steps

To get all necessary software for running `scprocess` load the Miniforge3 module by running

```bash
ml Miniforge3
```

## Installation

1. Clone the repository:

    ```bash
    git clone https://github.com/marusakod/scprocess.git
    ```

2.  Create a Conda environment:
    
    Navigate into the `scprocess` directory and create a Conda environment:
    
    ```bash
    conda env create -n scprocess -f envs/scprocess_slurm.yaml
    ```

3.  Add `scprocess` to your PATH. 
    
    Open your `~/.bashrc` file and add the following line:

    ```bash
    export PATH=/PATH/TO/YOUR/FOLDER/scprocess:${PATH}
    ```
    
    Verify the installation by:

    * Reloading your `~/.bashrc` file: 
      ```bash
      source ~/.bashrc
      ```
    
    * Activating the `scprocess` Conda environment and checking for help messages:
      ```bash
      conda activate scprocess
      scprocess
      ```

## `scprocess` data directory setup

`scprocess` requires a dedicated directory to store all necessary data, such as reference genomes.

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
      profile: slurm_shpc
    ref_txomes:
      tenx:
      - name: human_2024 
      - name: mouse_2024 
    ```

    This will ask the setup process to download and prepare the most recent pre-built [human](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads:~:text=Human%20reference%20(GRCh38)%20%2D%202024%2DA) and [mouse](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads:~:text=Mouse%20reference%20(GRCm39)%20%2D%202024%2DA) reference transcriptomes from 10x Genomics.
    
    Optionally you can add parameters which will be used to create a template configuration file using the `scprocess newproj -c` command for example:
    
    ```yaml
    user:
      profile: slurm_shpc
      your_name: Testy McUser
      affiliation: F. Hoffmann-La Roche Ltd.
    arvados:
      arv_instance: arkau
    ref_txomes:
      tenx:
      - name: human_2024 
      - name: mouse_2024 
    ```
    
    
    For more information on how to structure the _scprocess_setup.yaml_ see the [`Reference`](https://marusakod.github.io/scprocess/reference/) section.

4. Finish setting up the `scprocess` data directory:

    To download all required data and index reference transcriptomes use the `scprocess setup` command. The first time you run `scprocess setup` you need to specify a `-c`/`--rangerurl` flag and provide a valid download link for Cell Ranger (v9.0.0 or higher) available on the [10x Genomics Cell Ranger download & installation page](https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions) : 
    
    ```bash
    scprocess setup -c "https://cf.10xgenomics.com/releases/cell-exp/cellranger-10.0.0.tar.gz..." 
    ```
    
    Note that {{sc}} only requires barcode whitelists from Cell Ranger, therefore the full Cell Ranger installation will not be retained after the setup process.

    Once the inital setup is complete, you do not need to provide the Cell Ranger link again i.e if you modify the _scprocess_setup.yaml_ file to add additional reference genomes, simply run `scprocess setup`.

# Using `scprocess`

If you are using scprocess for the first time, we recommend working through the [Tutorials](https://marusakod.github.io/scprocess/tutorials/) and [Usage](https://marusakod.github.io/scprocess/usage/) sections of the documentation. For a detailed explanation of all parameters, please refer to the [Reference](https://marusakod.github.io/scprocess/reference/) page.

# sHPC-specific usage notes

## Accessing raw data from Arvados

`scprocess` supports raw data stored either in a local directory (defined via the `fastq_dir` parameter) or hosted on **Arvados**. To use Arvados directly, provide a list of collection UUIDs to the `arv_uuids` parameter as well as the name of the Arvados instance to the `arv_instance` parameter in your configuration file e.g:

```yaml
project: 
  arv_uuids: ["arkau-qr8st-1a2b3c4d5e6f7g8", "arkau-9v0wx-h9i8j7k6l5m4n3o", "arkau-z2y3x-p0q1r2s3t4u5v6w"]
  arv_instance: arkau
```
 
## Running Cellbender

The default ambient correction method in `scprocess` is DexontX which doesn't require any additional software. However, if you choose to use CellBender, Apptainer is required.

Because Apptainer is not available on login nodes, you must execute `scprocess` from an interactive session or via a batch job: 

1. Example 1: run `scprocess` from an interactive session:

  ```bash
  # start a tmux session. For more information about tmux visit https://github.com/tmux/tmux/wiki
  tmux new -s scprocess

  # start interactive session
  srun --pty --qos=interactive --partition=interactive_cpu -t 0-24:00 bash -l

  # set up environment
  ml Miniforge3
  conda activate scprocess

  # run scprocess
  scprocess run /path/to/config/config-project.yaml
  ```

2. example 2: submit batch job to run `scprocess`

  ```bash
  #!/bin/bash 
  #SBATCH -n 1                               
  #SBATCH --cpus-per-task=1                 
  #SBATCH -t 0-03:00                       

  ml purge
  ml Miniforge3
  conda activate scprocess

  scprocess run /path/to/config/config-project.yaml
  ```

