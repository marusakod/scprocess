
# How to set up `scprocess` on sHPC

## First steps

To get all necessary software for running `scprocess` load the Miniforge3 module by running

```bash
source /apps/rocs/init.sh 2020.08
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
      scprocess -h
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
      profile: slurm_shpc # scprocess will run with a job scheduler
      your_name:    Testy McUser # edit this
      affiliation:  Unemployed # edit this
    genomes:
      tenx:
      - name: human_2024 
      - name: mouse_2024 
    ```

    This will ask the setup process to download and prepare the most recent pre-built [human](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads:~:text=Human%20reference%20(GRCh38)%20%2D%202024%2DA) and [mouse](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads:~:text=Mouse%20reference%20(GRCm39)%20%2D%202024%2DA) genomes from 10x Genomics. For more information on how to structure the _scprocess_setup.yaml_ see the [`Reference`](https://macnairw.pages.roche.com/scprocess/reference)section.
    

4. Finish setting up the `scprocess` data directory with:
    ```bash
    scprocess setup 
    ```

## Advanced setup and parameters

### Accessing raw data from arvados

### Accessing raw data from Arvados

`scprocess` supports raw data stored either in a local directory (defined via the `fastq_dir` parameter) or hosted on **Arvados**. To use Arvados directly, provide a list of collection UUIDs to the `arv_uuids` parameter in your configuration file e.g:

```yaml
project: 
  arv_uuids: ["arkau-qr8st-1a2b3c4d5e6f7g8", "arkau-9v0wx-h9i8j7k6l5m4n3o", "arkau-z2y3x-p0q1r2s3t4u5v6w"]
```
    
To enable Arvados integration, you must [configure your Arvados API credentials](https://pdc.pages.roche.com/docs/service/pdc-home/get-started/arv-setup.html) and initialize the environment before execution:
    
```bash
ml purge
module load arvados
arv-env arkau
module load Miniforge3
conda activate scprocess
unset conda 
```

For convenience, you can add the following helper function to your `~/.bashrc` to automate this setup: 

```bash
function arv-setup() {
  # 1. Load modules
  module load arvados || { echo "Failed to load arvados"; return 1; }

  # 2. Set environment
  arv-env arkau

  # 3. Load Miniforge
  module load Miniforge3 || { echo "Failed to load Miniforge3"; return 1; }

  # 4. Activate Conda environment
  # Using 'source' ensures the shell picks up the environment correctly
  conda activate scprocess

  unset conda

  echo "Environment 'scprocess' is now active with Arvados tools."
}
```

Once added, simply run `arv-setup` in your terminal to prepare your session.

### Running Cellbender

The default ambient correction method in `scprocess` is DexontX which doesn't require any additional software. However, if you choose to use CellBender, Apptainer is required.

Because Apptainer is not available on login nodes, you must execute `scprocess` from an interactive session or via a batch job: 

1. Example 1: run `scprocess` from an interactive session:

  ```bash
  # start a tmux session. For more information about tmux visit https://github.com/tmux/tmux/wiki
  tmux new -s scprocess

  # start interactive session
  srun --pty --qos=interactive --partition=interactive_cpu -t 0-24:00 bash -l

  # set up environment
  arv-setup

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

