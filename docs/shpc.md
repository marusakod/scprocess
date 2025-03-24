# How to set up {{ software_name }} on the `sHPC`

To get {{ software_name }} to work, we need to connect it up to the functionality available on the `sHPC`.

1. Clone the {{ software_name }} repo to your home directory.

2. Change to the main-shpc branch:

```
git checkout main-shpc
git branch
#   dev-hto
#   main
# * main-shpc
```

You should be able to see a file _lsf.yaml_ in the top level of the {{ software_name }} directory.

3. Add some things to your `.bashrc`:

```
alias scprocess='export ROCS_ARCH=sandybridge; source /apps/rocs/init.sh; ml snakemake-lsf/1.0.7-foss-2020a-Python-3.11.3-snakemake-8.23.0; scprocess'
alias scsetup='export ROCS_ARCH=sandybridge; source /apps/rocs/init.sh; ml snakemake-lsf/1.0.7-foss-2020a-Python-3.11.3-snakemake-8.23.0; scsetup'
```
