# Plan, fit and evaluate an XGBoost cell-type classifier.

if TRAIN_XGB_PARAMS is not None:

  # Hold out whole groups and select cells before loading expression matrices.
  rule train_xgboost_prepare:
    input:
      annots_f    = TRAIN_XGB_PARAMS['annots_f'],
      integration_f = f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
      h5ads_yaml  = f'{int_dir}/h5ads_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml',
      label_map   = [TRAIN_XGB_PARAMS['label_map_f']] if TRAIN_XGB_PARAMS.get('label_map_f') else [],
      planner     = str(scprocess_dir / 'scripts/train_xgboost.py'),
    output:
      training_metadata_f     = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_training_metadata.csv.gz',
      coverage_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_split_coverage.csv',
      counts_f   = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_label_counts.csv',
      metadata_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_training_metadata.csv.gz.json',
    params:
      excluded_labels                  = TRAIN_XGB_PARAMS['excluded_labels'],
      label_map_f                      = TRAIN_XGB_PARAMS.get('label_map_f') or '',
      cluster_majority_vote_resolution = TRAIN_XGB_PARAMS.get('cluster_majority_vote_resolution', ''),
      cluster_majority_vote_min_purity = TRAIN_XGB_PARAMS['cluster_majority_vote_min_purity'],
      test_group_var                   = TRAIN_XGB_PARAMS.get('test_group_var') or '',
      test_groups                      = TRAIN_XGB_PARAMS.get('test_groups', []),
      validation_group_var             = TRAIN_XGB_PARAMS['validation_group_var'],
      validation_fraction              = TRAIN_XGB_PARAMS['validation_fraction'],
      validation_variable_combinations = TRAIN_XGB_PARAMS['validation_variable_combinations'],
      sampling_variable_combinations   = TRAIN_XGB_PARAMS['sampling_variable_combinations'],
      train_cells_per_label             = TRAIN_XGB_PARAMS['train_cells_per_label'],
      validation_cells_per_label        = TRAIN_XGB_PARAMS['validation_cells_per_label'],
      min_cells_per_label               = TRAIN_XGB_PARAMS['min_cells_per_label'],
      seed                             = TRAIN_XGB_PARAMS['seed'],
      batch_var                        = BATCH_VAR,
      is_join                          = 0,
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb  = config['resources']['gb_train_xgboost_prepare'] * 1024,
      runtime = config['resources']['mins_train_xgboost_prepare']
    log: f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_training_metadata.csv.gz.log'
    conda: '../envs/label_celltypes.yaml'
    shell: """
      exec &>> {log}
      python3 {scprocess_dir}/scripts/train_xgboost.py prepare \
        --annots_f                          {input.annots_f:q} \
        --integration_f                       {input.integration_f:q} \
        --h5ads_yaml                        {input.h5ads_yaml:q} \
        --batch_var                         {params.batch_var:q} \
        --is_join                           {params.is_join} \
        --excluded_labels                   {params.excluded_labels:q} \
        --label_map_f={params.label_map_f:q} \
        --cluster_majority_vote_resolution={params.cluster_majority_vote_resolution:q} \
        --cluster_majority_vote_min_purity  {params.cluster_majority_vote_min_purity:q} \
        --test_group_var={params.test_group_var:q} \
        --test_groups                       {params.test_groups:q} \
        --validation_group_var              {params.validation_group_var:q} \
        --validation_fraction               {params.validation_fraction:q} \
        --validation_variable_combinations  {params.validation_variable_combinations:q} \
        --sampling_variable_combinations    {params.sampling_variable_combinations:q} \
        --train_cells_per_label              {params.train_cells_per_label:q} \
        --validation_cells_per_label         {params.validation_cells_per_label:q} \
        --min_cells_per_label                {params.min_cells_per_label:q} \
        --seed                              {params.seed:q} \
        --training_metadata_f                            {output.training_metadata_f:q} \
        --coverage_f                        {output.coverage_f:q} \
        --counts_f                          {output.counts_f:q}
      """


  # Fit both passes and export the selected genes and validation predictions.
  rule train_xgboost_train:
    input:
      training_metadata_f     = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_training_metadata.csv.gz',
      h5ads_yaml = f'{int_dir}/h5ads_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml',
      trainer    = str(scprocess_dir / "scripts/train_xgboost.py"),
    output:
      model_f      = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_xgboost_model.json',
      cls_f        = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_allowed_cls.csv',
      genes_f      = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_selected_genes.txt',
      imp_f        = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_gene_importance.csv',
      preds_f      = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_predictions.csv.gz',
      pb_f         = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_pseudobulk.h5ad',
      model_card_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_model_card.json',
    params:
      output_dir             = xgb_dir,
      xgboost_tag            = TRAIN_XGB_PARAMS['xgboost_tag'],
      use_gpu                = int(TRAIN_XGB_PARAMS['use_gpu']),
      gene_exclude_re        = TRAIN_XGB_PARAMS.get('gene_exclude_re') or '',
      min_cells_expressed    = TRAIN_XGB_PARAMS['min_cells_expressed'],
      seed                   = TRAIN_XGB_PARAMS['seed'],
      pass1_subsample        = TRAIN_XGB_PARAMS['pass1_subsample'],
      pass1_colsample_bytree = TRAIN_XGB_PARAMS['pass1_colsample_bytree'],
      pass1_learning_rate    = TRAIN_XGB_PARAMS['pass1_learning_rate'],
      pass1_nrounds          = TRAIN_XGB_PARAMS['pass1_nrounds'],
      pass1_early_stopping   = TRAIN_XGB_PARAMS['pass1_early_stopping'],
      gain_threshold         = TRAIN_XGB_PARAMS['gain_threshold'],
      min_genes              = TRAIN_XGB_PARAMS['min_genes'],
      max_genes              = TRAIN_XGB_PARAMS['max_genes'],
      pass2_colsample_bytree = TRAIN_XGB_PARAMS['pass2_colsample_bytree'],
      pass2_learning_rate    = TRAIN_XGB_PARAMS['pass2_learning_rate'],
      pass2_nrounds          = TRAIN_XGB_PARAMS['pass2_nrounds'],
      pass2_early_stopping   = TRAIN_XGB_PARAMS['pass2_early_stopping'],
      label_map_f            = TRAIN_XGB_PARAMS.get('label_map_f') or '',
    threads: 8
    retries: config['resources']['retries']
    resources:
      mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'train_xgboost_train', 'memory', attempt),
      runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'train_xgboost_train', 'time', attempt)
    log: f'{logs_dir}/train_xgboost/train_xgboost_{TRAIN_XGB_PARAMS["xgboost_tag"]}_{DATE_STAMP}.log'
    benchmark: f'{benchmark_dir}/train_xgboost/train_xgboost_{TRAIN_XGB_PARAMS["xgboost_tag"]}_{DATE_STAMP}.benchmark.txt'
    conda:
      '../envs/train_xgboost_gpu.yaml' if TRAIN_XGB_PARAMS['use_gpu'] else '../envs/label_celltypes.yaml'
    shell: """
      exec &>> {log}
      python3 {scprocess_dir}/scripts/train_xgboost.py train \
        --training_metadata_f                            {input.training_metadata_f:q} \
        --h5ads_yaml                        {input.h5ads_yaml:q} \
        --xgboost_tag                       {params.xgboost_tag:q} \
        --use_gpu                           {params.use_gpu:q} \
        --gene_exclude_re={params.gene_exclude_re:q} \
        --min_cells_expressed               {params.min_cells_expressed:q} \
        --seed                              {params.seed:q} \
        --pass1_subsample                   {params.pass1_subsample:q} \
        --pass1_colsample_bytree            {params.pass1_colsample_bytree:q} \
        --pass1_learning_rate               {params.pass1_learning_rate:q} \
        --pass1_nrounds                     {params.pass1_nrounds:q} \
        --pass1_early_stopping              {params.pass1_early_stopping:q} \
        --gain_threshold                    {params.gain_threshold:q} \
        --min_genes                         {params.min_genes:q} \
        --max_genes                         {params.max_genes:q} \
        --pass2_colsample_bytree            {params.pass2_colsample_bytree:q} \
        --pass2_learning_rate               {params.pass2_learning_rate:q} \
        --pass2_nrounds                     {params.pass2_nrounds:q} \
        --pass2_early_stopping              {params.pass2_early_stopping:q} \
        --label_map_f={params.label_map_f:q} \
        --output_dir                        {params.output_dir:q} \
        --n_cores                           {threads}
      """


  # Apply the fitted model to every cell in each batch.
  rule train_xgboost_predict:
    input:
      adata_f = f'{int_dir}/anndata_cells_clean_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5ad',
      model_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_xgboost_model.json',
      cls_f   = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_allowed_cls.csv',
      genes_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_selected_genes.txt',
    output:
      pred_f = temp(f'{xgb_dir}/tmp_fulldata_predict_{TRAIN_XGB_PARAMS["xgboost_tag"]}_{{batch}}.csv.gz')
    params:
      ref_tag   = TRAIN_XGB_PARAMS['xgboost_tag'],
      batch_var = TRAIN_XGB_PARAMS.get('batch_var', BATCH_VAR)
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'train_xgboost_predict', 'memory', attempt),
      runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'train_xgboost_predict', 'time', attempt)
    log: f'{logs_dir}/train_xgboost/train_xgboost_predict_{TRAIN_XGB_PARAMS["xgboost_tag"]}_{{batch}}_{DATE_STAMP}.log'
    benchmark: f'{benchmark_dir}/train_xgboost/train_xgboost_predict_{TRAIN_XGB_PARAMS["xgboost_tag"]}_{{batch}}_{DATE_STAMP}.benchmark.txt'
    conda:
      '../envs/label_celltypes.yaml'
    shell: """
      exec &>> {log}
      python3 {scprocess_dir}/scripts/label_celltypes.py xgboost_one_batch \
        {wildcards.batch} {params.batch_var} {params.ref_tag} \
        --adata_f                           {input.adata_f} \
        --model_f                           {input.model_f} \
        --cls_f                             {input.cls_f} \
        --genes_f                           {input.genes_f} \
        --pred_f                            {output.pred_f}
      """


  # Apply display specifications to saved predictions without retraining.
  rule train_xgboost_specs:
    input:
      cls_f             = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_allowed_cls.csv',
      counts_f          = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_label_counts.csv',
      predictions_f     = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_fulldata_predictions.csv.gz',
      specifications    = get_xgboost_spec_inputs(TRAIN_XGB_PARAMS),
      script_f          = str(scprocess_dir / "scripts/train_xgboost.py"),
    output:
      label_specs_f                 = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_label_specs.csv',
      coarse_label_specs_f          = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_coarse_label_specs.csv',
      coarse_evaluation_f           = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_coarse_evaluation.csv',
      coarse_per_label_evaluation_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_coarse_per_label_evaluation.csv',
    params:
      xgboost_tag          = TRAIN_XGB_PARAMS['xgboost_tag'],
      label_specs_f        = TRAIN_XGB_PARAMS.get('label_specs_f') or '',
      coarse_label_specs_f = TRAIN_XGB_PARAMS.get('coarse_label_specs_f') or '',
      label_map_f          = TRAIN_XGB_PARAMS.get('label_map_f') or '',
      output_dir           = str(pathlib.Path(f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}').parent),
    threads: 1
    resources:
      mem_mb  = 2048,
      runtime = 5
    log: f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_specs.log'
    conda: '../envs/label_celltypes.yaml'
    shell: """
      exec &>> {log}
      python3 {scprocess_dir}/scripts/train_xgboost.py specs \
        --cls_f                             {input.cls_f:q} \
        --counts_f                          {input.counts_f:q} \
        --predictions_f                     {input.predictions_f:q} \
        --xgboost_tag                       {params.xgboost_tag:q} \
        --label_specs_f={params.label_specs_f:q} \
        --coarse_label_specs_f={params.coarse_label_specs_f:q} \
        --label_map_f={params.label_map_f:q} \
        --output_dir                        {params.output_dir:q}
      """


  # Combine full-data calls and evaluate each saved split.
  rule train_xgboost_combine:
    input:
      pred_fs = expand(
        f'{xgb_dir}/tmp_fulldata_predict_{TRAIN_XGB_PARAMS["xgboost_tag"]}_{{batch}}.csv.gz',
        batch = BATCHES),
      training_metadata_f      = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_training_metadata.csv.gz',
      integration_f = f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    output:
      fulldata_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_fulldata_predictions.csv.gz'
    params:
      label_map_f = TRAIN_XGB_PARAMS.get('label_map_f', ''),
    retries: config['resources']['retries']
    resources:
      mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'train_xgboost_combine', 'memory', attempt),
      runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'train_xgboost_combine', 'time', attempt)
    log: f'{logs_dir}/train_xgboost/train_xgboost_combine_{TRAIN_XGB_PARAMS["xgboost_tag"]}_{DATE_STAMP}.log'
    benchmark: f'{benchmark_dir}/train_xgboost/train_xgboost_combine_{TRAIN_XGB_PARAMS["xgboost_tag"]}_{DATE_STAMP}.benchmark.txt'
    conda:
      '../envs/label_celltypes.yaml'
    shell: """
      exec &>> {log}
      python3 {scprocess_dir}/scripts/train_xgboost.py combine \
        {input.pred_fs:q} --training_metadata_f {input.training_metadata_f:q} --integration_f {input.integration_f:q} \
        --output_f                          {output.fulldata_f:q}
      """
