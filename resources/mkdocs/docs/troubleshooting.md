# Troubleshooting

## ~~Frequently~~ Possibly Asked Questions

??? question "I am having issues with setting up the {{sc}} conda environment"

    Try deleting your _~/.condarc_ file.

??? question "I get an error about `snakemake cannot lock directory` when I try to call {{scrun}}"

    To avoid nasty surprises, {{sc}} locks several directories before running. If you cancel a {{sc}} run (or it otherwise gets stopped unexpectedly), this locking can remain. To fix this, you can run `scprocess run config-my_project.yaml --unlock`. This should give you a message about `snakemake unlocked the directory`, and then you run `scprocess run config-my_project.yaml` as normal.

??? question "One of my rules times out"

    You can change the default parameters for a rule by editing the relevant profile, e.g. _profile/slurm_default/config.yaml_ in your {{sc}} directory. For example, to increase the time allowed for the run_mapping and run_mapping_hto rules on our system, I added these lines:
    ```
    set-resources:
      run_mapping:
        qos: '1d'
        runtime: 1440
      run_mapping_hto:
        qos: '1d'
        runtime: 1440
    ```
    To increase the number of threads available to these rules, I added these lines:
    ```
    set-resources:
      run_mapping: 32
      run_mapping_hto: 32
    ```


## Bug reporting

If you encounter any problems when running {{sc}} please consider opening an issue on our [GitHub](https://github.com/marusakod/scprocess).
