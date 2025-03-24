# How to set up {{ software_name }} on the `sHPC`

To get {{ software_name }} to work, we need download it and link it to the functionality available on the `sHPC`.

1. Clone the {{ software_name }} repo to your home directory.

```bash
cd ~/packages/ # or wherever you keep your packages
git clone https://code.roche.com/macnairw/scprocess
```

2. Change to the `main-shpc` branch:

```bash
# change to relevant branch
git checkout main-shpc
# check it worked (the important thing is that it says * main-shpc; the other names might be different)
git branch
#   dev-hto
#   main
# * main-shpc
```

You should be able to see a file _lsf.yaml_ in the top level of the {{ software_name }} directory:

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

3. Add some things to your `~/.bashrc`:

```bash
# add scprocess to path
echo "export PATH=~/packages/scprocess:${PATH}" >> ~/.bashrc
# add some sHPC-specific things
echo "alias scprocess='export ROCS_ARCH=sandybridge; source /apps/rocs/init.sh; ml snakemake-lsf/1.0.7-foss-2020a-Python-3.11.3-snakemake-8.23.0; scprocess'" >> ~/.bashrc
echo "alias scsetup='export ROCS_ARCH=sandybridge; source /apps/rocs/init.sh; ml snakemake-lsf/1.0.7-foss-2020a-Python-3.11.3-snakemake-8.23.0; scsetup'" >> ~/.bashrc
# this code adds some extra lines to the end of your .bashrc file. feel free to put them somewhere more tidy!
```

Check that this worked:

```bash
# reload the .bashrc file
source ~/.bashrc

# check that scprocess works
scprocess -h
```

4. You're now ready go to [Getting started](setup.md).
