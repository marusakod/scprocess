# Create every conda environment used by scprocess rules, without running them.
# Called by `scprocess setup --create-envs` with --conda-create-envs-only, so the
# placeholder outputs below are never written. Snakemake identifies each
# environment by its file contents and conda prefix, so later project and join
# runs using the same prefix reuse these environments.

# environment files referenced by `conda:` directives in rules/*.smk;
# tests/test_create_envs.py checks that this list stays complete
RULE_ENV_NAMES = [
  'alevin_fry',
  'ambientr',
  'bpcells_pca',
  'hvgs',
  'integration',
  'label_celltypes',
  'py_env',
  'rlibs',
  'rlibs_bpcells',
  'scprocess_local',
  'shiny',
  'train_xgboost_gpu',
  'tricycle'
]


rule create_envs:
  input:
    [f'.scprocess_create_env_{env_name}' for env_name in RULE_ENV_NAMES]
  default_target: True


# one placeholder rule per environment, e.g. create_env_rlibs uses envs/rlibs.yaml
for env_name in RULE_ENV_NAMES:
  rule:
    name:
      f'create_env_{env_name}'
    output:
      f'.scprocess_create_env_{env_name}'
    conda:
      f'../envs/{env_name}.yaml'
    shell:
      'touch {output}'
