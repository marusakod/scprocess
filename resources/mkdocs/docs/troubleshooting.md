# Troubleshooting

## ~~Frequently~~ Possibly Asked Questions

??? question "I changed some parameters, and I would like to rerun {{sc}}"

    One slightly hacky way of doing this is to delete the first output file that would be affected by the parameter change. For example, if you ran {{scrun}} and it completed, and you would like to update the results after changing a QC parameter, then you could delete the file _output/my_project_qc/qc_dt_all_samples_my_project_2025-01-01.txt.gz_. Then you can run `scprocess run config-my_project.yaml -n` and you should see only relevant steps proposed.

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
