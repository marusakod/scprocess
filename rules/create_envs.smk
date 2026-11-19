# Create every conda environment used by scprocess rules, without running them.
# Called by `scprocess setup --create-envs`

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
