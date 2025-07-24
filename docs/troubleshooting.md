# Troubleshooting

## ~~Frequently~~ Possibly Asked Questions

??? I changed some parameters, and I would like to rerun `scprocess`
  One slightly hacky way of doing this is to delete the first output file that would be affected by the parameter change. For example, if you ran `scprocess run` and it completed, and you would like to update the results after changing a QC parameter, then you could delete the file _output/my_project_qc/qc_dt_all_samples_my_project_2025-01-01.txt.gz_. Then you can run `scprocess run config-my_project.yaml -n` and you should see only relevant steps proposed.

??? I get an error about `snakemake cannot lock directory` when I try to call `scprocess run`
  To avoid nasty surprises, `snakemake` locks several directories before running. If you cancel a `snakemake` run (or it otherwise gets stopped unexpectedly), this locking can remain. To fix this, you can run `scprocess run config-my_project.yaml --unlock`. This should give you a message about `snakemake unlocked the directory`, and then you run `scprocess run config-my_project.yaml` as normal.


## Bug reporting

If you encounter any problems when running {{ software_name }} please consider opening an issue on our [GitHub](https://github.com/marusakod/scprocess_test).
