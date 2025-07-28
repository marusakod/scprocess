# Troubleshooting

## ~~Frequently~~ Possibly Asked Questions

??? question "I looked at the knee plots and the knees and shins look wrong for some of them"

    We're sorry to hear that. Often this is a sign that those samples are bad quality :( But at least we make it relative straightforward to handle this, by specifying the correct knees and shins by hand. 

    You need to specify them in a file containing custom parameters that you need to create, e.g. _data/metadata/custom_params_my_project_2025-01-01.yaml_. What we recommend doing is to use {{scknee}}, e.g. by calling `scprocess plotknee -s id_of_bad_sample -c config-my_project.yaml`. This will create an html version of the knee plot, in mapping directory for that sample (e.g. in _output/my_project_mapping/af_id_of_bad_sample_). If you open this and hover your mouse over the correct points of the curve, to find the right number of UMIs for `knee1`, `shin1`, `knee2` and `shin2`. 

    You can put these values into your custom parameters yaml file, which might look like this:
    ```yaml
    id_of_bad_sample:
      ambient:
        knee1: 2000
        shin1: 1000
        knee2: 150
        shin2: 50
    ```
    Finally, you will need to tell {{sc}} to use this file, by pointing to it in your config yaml file. If _config-myproject.yaml_ is your file, you'll need to add a line like this:
    ```yaml
    custom_sample_params: data/metadata/custom_params_my_project_2025-01-01.yaml
    ```
    It doesn't matter where you put it, but a nice place is just above the `ambient` block.


??? question "I changed some parameters, and I would like to rerun {{sc}}"

    One slightly hacky way of doing this is to delete the first output file that would be affected by the parameter change. For example, if you ran {{scrun}} and it completed, and you would like to update the results after changing a QC parameter, then you could delete the file _output/my_project_qc/qc_dt_all_samples_my_project_2025-01-01.txt.gz_. Then you can run `scprocess run config-my_project.yaml -n` and you should see only relevant steps proposed.

??? question "I get an error about `snakemake cannot lock directory` when I try to call {{scrun}}"

    To avoid nasty surprises, {{sc}} locks several directories before running. If you cancel a {{sc}} run (or it otherwise gets stopped unexpectedly), this locking can remain. To fix this, you can run `scprocess run config-my_project.yaml --unlock`. This should give you a message about `snakemake unlocked the directory`, and then you run `scprocess run config-my_project.yaml` as normal.


## Bug reporting

If you encounter any problems when running {{sc}} please consider opening an issue on our [GitHub](https://github.com/marusakod/scprocess_test).
