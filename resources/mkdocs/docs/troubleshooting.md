# Troubleshooting

## ~~Frequently~~ Possibly Asked Questions

??? question "I get an error about `snakemake cannot lock directory` when I try to call {{scrun}}"

    To avoid nasty surprises, {{sc}} locks several directories before running. If you cancel a {{sc}} run (or it otherwise gets stopped unexpectedly), this locking can remain. To fix this, you can run `scprocess run config-my_project.yaml --unlock`. This should give you a message about `snakemake unlocked the directory`, and then you run `scprocess run config-my_project.yaml` as normal.

??? question "One of my rules times out/fails after reaching memory limit"

    When a `Snakemake` rule fails, the error message typically specifies the exact rule name, for example: `Error in rule run_mapping`. To resolve resource-related errors you can override the default resource limits in your configuration file. For instance, to allocate 16 GB of memory and 60 minutes of runtime to the `run_mapping` rule, add the following to your config:

    ```yaml
    resources:
      gb_run_mapping: 16 # memory in GB
      mins_run_mapping: 60 # runtime in minutes
    ```
    For more information please see the [Reference](reference.md#resources) section.

## Bug reporting

If you encounter any problems when running {{sc}} please consider opening an issue on our [GitHub](https://github.com/marusakod/scprocess).


