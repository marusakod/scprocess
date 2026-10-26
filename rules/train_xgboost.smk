import json
if TRAIN_XGB_PARAMS is not None:

  rule train_xgboost_plan:
    input:
      annots_f = TRAIN_XGB_PARAMS['annots_f'],
      cluster_csv = f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
      h5ads_yaml = f'{int_dir}/h5ads_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml',
      label_map = [TRAIN_XGB_PARAMS['label_map_f']] if TRAIN_XGB_PARAMS.get('label_map_f') else [],
      planner = str(scprocess_dir / 'scripts/train_xgboost_plan.py'),
    output:
      plan_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_training_plan.csv.gz',
      coverage_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_split_coverage.csv',
      counts_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_label_counts.csv',
      metadata_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_training_plan.csv.gz' + '.json',
    params:
      config_json = json.dumps(get_xgboost_planning_config(TRAIN_XGB_PARAMS)),
      batch_var = BATCH_VAR,
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb = config['resources'].get('gb_train_xgboost_plan', 8) * 1024,
      runtime = config['resources'].get('mins_train_xgboost_plan', 10)
    log: f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_training_plan.csv.gz' + '.log'
    conda: '../envs/label_celltypes.yaml'
    shell: """
      exec &>> {log}
      python3 {scprocess_dir}/scripts/train_xgboost.py plan \
        --annots_f {input.annots_f:q} --cluster_csv {input.cluster_csv:q} \
        --h5ads_yaml {input.h5ads_yaml:q} --batch_var {params.batch_var:q} \
        --config_json {params.config_json:q} --plan_f {output.plan_f:q} \
        --coverage_f {output.coverage_f:q} --counts_f {output.counts_f:q}
      """


  rule train_xgboost_train:
    input:
      plan_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_training_plan.csv.gz',
      h5ads_yaml = f'{int_dir}/h5ads_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml',
      trainer = str(scprocess_dir / "scripts/train_xgboost_workflow.py"),
    output:
      model_f   = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_xgboost_model.json',
      cls_f     = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_allowed_cls.csv',
      genes_f   = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_selected_genes.txt',
      imp_f     = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_gene_importance.csv',
      preds_f   = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_predictions.csv.gz',
      pb_f      = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_pseudobulk.h5ad',
      model_card_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_model_card.json',
    params:
      output_dir = xgb_dir,
      config_json = json.dumps(TRAIN_XGB_PARAMS),
    threads: 8
    retries: config['resources']['retries']
    resources:
      mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'train_xgboost_train', 'memory', attempt),
      runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'train_xgboost_train', 'time', attempt)
    log: f'{logs_dir}/train_xgboost/train_xgboost_{TRAIN_XGB_PARAMS["xgboost_tag"]}_{DATE_STAMP}.log'
    benchmark: f'{benchmark_dir}/train_xgboost/train_xgboost_{TRAIN_XGB_PARAMS["xgboost_tag"]}_{DATE_STAMP}.benchmark.txt'
    conda:
      '../envs/train_xgboost_gpu.yaml' if TRAIN_XGB_PARAMS.get('use_gpu', False) else '../envs/label_celltypes.yaml'
    shell: """
      exec &>> {log}
      python3 {scprocess_dir}/scripts/train_xgboost.py train \
        --plan_f {input.plan_f:q} --h5ads_yaml {input.h5ads_yaml:q} \
        --config_json {params.config_json:q} --output_dir {params.output_dir:q} \
        --n_cores {threads}
      """


  rule train_xgboost_predict:
    input:
      adata_f = f'{int_dir}/anndata_cells_clean_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5ad',
      model_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_xgboost_model.json',
      cls_f   = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_allowed_cls.csv',
      genes_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_selected_genes.txt',
    output:
      pred_f  = temp(f'{xgb_dir}/tmp_fulldata_predict_{TRAIN_XGB_PARAMS["xgboost_tag"]}_{{batch}}.csv.gz')
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
        --adata_f   {input.adata_f} \
        --model_f   {input.model_f} \
        --cls_f     {input.cls_f} \
        --genes_f   {input.genes_f} \
        --pred_f    {output.pred_f}
      """


  rule train_xgboost_aggregate:
    input:
      pred_fs           = expand(
        f'{xgb_dir}/tmp_fulldata_predict_{TRAIN_XGB_PARAMS["xgboost_tag"]}_{{batch}}.csv.gz',
        batch=BATCHES),


      plan_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_training_plan.csv.gz',
      cluster_csv = f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    output:
      fulldata_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["xgboost_tag"]}_fulldata_predictions.csv.gz'
    params:
      label_map_f = TRAIN_XGB_PARAMS.get('label_map_f', ''),
    retries: config['resources']['retries']
    resources:
      mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'train_xgboost_aggregate', 'memory', attempt),
      runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'train_xgboost_aggregate', 'time', attempt)
    log: f'{logs_dir}/train_xgboost/train_xgboost_aggregate_{TRAIN_XGB_PARAMS["xgboost_tag"]}_{DATE_STAMP}.log'
    benchmark: f'{benchmark_dir}/train_xgboost/train_xgboost_aggregate_{TRAIN_XGB_PARAMS["xgboost_tag"]}_{DATE_STAMP}.benchmark.txt'
    conda:
      '../envs/label_celltypes.yaml'
    shell: """
      exec &>> {log}
      python3 {scprocess_dir}/scripts/train_xgboost.py aggregate_fulldata \
        {input.pred_fs:q} --plan_f {input.plan_f:q} --cluster_csv {input.cluster_csv:q} \
        --output_f {output.fulldata_f:q}
      """
