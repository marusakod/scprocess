if TRAIN_XGB_PARAMS is not None:

  rule train_xgboost_train:
    input:
      annots_f    = TRAIN_XGB_PARAMS['annots_f'],
      cluster_csv = f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
      h5ads_yaml  = f'{int_dir}/h5ads_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml',
    output:
      model_f   = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_xgboost_model.json',
      cls_f     = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_allowed_cls.csv',
      genes_f   = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_selected_genes.txt',
      imp_f     = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_gene_importance.csv',
      preds_f   = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_predictions.csv.gz',
      pb_f      = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_pseudobulk.h5ad',
      label_counts_f  = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_label_counts.csv',
    params:
      ref_tag              = TRAIN_XGB_PARAMS['ref_tag'],
      output_dir           = xgb_dir,
      batch_var            = TRAIN_XGB_PARAMS.get('batch_var', BATCH_VAR),
      int_res_ls           = TRAIN_XGB_PARAMS.get('int_res_ls', config.get('integration', {}).get('int_res_ls', [0.1, 0.2, 0.5, 1, 2])),
      label_map_f          = TRAIN_XGB_PARAMS.get('label_map_f', ''),
      refine_labels        = TRAIN_XGB_PARAMS.get('refine_labels', True),
      purity_threshold     = TRAIN_XGB_PARAMS.get('purity_threshold', 0.65),
      n_cells_per_type     = TRAIN_XGB_PARAMS.get('n_cells_per_type', 1000),
      min_cells_per_type   = TRAIN_XGB_PARAMS.get('min_cells_per_type', 20),
      min_cells_expressed  = TRAIN_XGB_PARAMS.get('min_cells_expressed', 10),
      gene_exclude_re      = TRAIN_XGB_PARAMS.get('gene_exclude_re', '(lincRNA|lncRNA|pseudogene|antisense)'),
      seed                 = TRAIN_XGB_PARAMS.get('seed', 42),
      use_gpu              = TRAIN_XGB_PARAMS.get('use_gpu', False),
      pass1_subsample      = TRAIN_XGB_PARAMS.get('pass1_subsample', 0.632),
      pass1_colsample_bytree = TRAIN_XGB_PARAMS.get('pass1_colsample_bytree', 0.1),
      pass1_learning_rate  = TRAIN_XGB_PARAMS.get('pass1_learning_rate', 0.1),
      pass1_nrounds        = TRAIN_XGB_PARAMS.get('pass1_nrounds', 300),
      pass1_early_stopping = TRAIN_XGB_PARAMS.get('pass1_early_stopping', 10),
      pass2_colsample_bytree = TRAIN_XGB_PARAMS.get('pass2_colsample_bytree', 0.5),
      pass2_learning_rate  = TRAIN_XGB_PARAMS.get('pass2_learning_rate', 0.05),
      pass2_nrounds        = TRAIN_XGB_PARAMS.get('pass2_nrounds', 500),
      pass2_early_stopping = TRAIN_XGB_PARAMS.get('pass2_early_stopping', 10),
      gain_threshold       = TRAIN_XGB_PARAMS.get('gain_threshold', 0.9),
      min_genes            = TRAIN_XGB_PARAMS.get('min_genes', 100),
      max_genes            = TRAIN_XGB_PARAMS.get('max_genes', 3000),
    threads: 8
    retries: config['resources']['retries']
    resources:
      mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'train_xgboost_train', 'memory', attempt),
      runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'train_xgboost_train', 'time', attempt)
    log: f'{logs_dir}/train_xgboost/train_xgboost_{TRAIN_XGB_PARAMS["ref_tag"]}_{DATE_STAMP}.log'
    benchmark: f'{benchmark_dir}/train_xgboost/train_xgboost_{TRAIN_XGB_PARAMS["ref_tag"]}_{DATE_STAMP}.benchmark.txt'
    conda:
      '../envs/label_celltypes.yaml'
    shell: """
      exec &>> {log}
      python3 {scprocess_dir}/scripts/train_xgboost.py train \
        --annots_f          {input.annots_f} \
        --cluster_csv       {input.cluster_csv} \
        --h5ads_yaml        {input.h5ads_yaml} \
        --ref_tag           {params.ref_tag} \
        --output_dir        {params.output_dir} \
        --batch_var         {params.batch_var} \
        --int_res_ls        "{params.int_res_ls}" \
        --refine_labels     {params.refine_labels} \
        --purity_threshold  {params.purity_threshold} \
        --n_cells_per_type  {params.n_cells_per_type} \
        --min_cells_per_type {params.min_cells_per_type} \
        --min_cells_expressed {params.min_cells_expressed} \
        --gene_exclude_re  "{params.gene_exclude_re}" \
        --seed              {params.seed} \
        --n_cores           {threads} \
        --pass1_subsample      {params.pass1_subsample} \
        --pass1_colsample_bytree {params.pass1_colsample_bytree} \
        --pass1_learning_rate  {params.pass1_learning_rate} \
        --pass1_nrounds        {params.pass1_nrounds} \
        --pass1_early_stopping {params.pass1_early_stopping} \
        --pass2_colsample_bytree {params.pass2_colsample_bytree} \
        --pass2_learning_rate  {params.pass2_learning_rate} \
        --pass2_nrounds        {params.pass2_nrounds} \
        --pass2_early_stopping {params.pass2_early_stopping} \
        --gain_threshold    {params.gain_threshold} \
        --min_genes         {params.min_genes} \
        --max_genes         {params.max_genes} \
        $( [ "{params.label_map_f}" != "" ] && echo "--label_map_f {params.label_map_f}" ) \
        $( [ "{params.use_gpu}" == "True" ] && echo "--use_gpu" )
      """


  rule train_xgboost_predict:
    input:
      adata_f = f'{int_dir}/anndata_cells_clean_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5ad',
      model_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_xgboost_model.json',
      cls_f   = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_allowed_cls.csv',
      genes_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_selected_genes.txt',
    output:
      pred_f  = temp(f'{xgb_dir}/tmp_fulldata_predict_{TRAIN_XGB_PARAMS["ref_tag"]}_{{batch}}.csv.gz')
    params:
      ref_tag   = TRAIN_XGB_PARAMS['ref_tag'],
      batch_var = TRAIN_XGB_PARAMS.get('batch_var', BATCH_VAR)
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'train_xgboost_predict', 'memory', attempt),
      runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'train_xgboost_predict', 'time', attempt)
    log: f'{logs_dir}/train_xgboost/train_xgboost_predict_{TRAIN_XGB_PARAMS["ref_tag"]}_{{batch}}_{DATE_STAMP}.log'
    benchmark: f'{benchmark_dir}/train_xgboost/train_xgboost_predict_{TRAIN_XGB_PARAMS["ref_tag"]}_{{batch}}_{DATE_STAMP}.benchmark.txt'
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
        f'{xgb_dir}/tmp_fulldata_predict_{TRAIN_XGB_PARAMS["ref_tag"]}_{{batch}}.csv.gz',
        batch=BATCHES),
      subsample_preds_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_predictions.csv.gz',
      annots_f          = TRAIN_XGB_PARAMS['annots_f'],
    output:
      fulldata_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_fulldata_predictions.csv.gz'
    params:
      label_map_f = TRAIN_XGB_PARAMS.get('label_map_f', ''),
    retries: config['resources']['retries']
    resources:
      mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'train_xgboost_aggregate', 'memory', attempt),
      runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'train_xgboost_aggregate', 'time', attempt)
    log: f'{logs_dir}/train_xgboost/train_xgboost_aggregate_{TRAIN_XGB_PARAMS["ref_tag"]}_{DATE_STAMP}.log'
    benchmark: f'{benchmark_dir}/train_xgboost/train_xgboost_aggregate_{TRAIN_XGB_PARAMS["ref_tag"]}_{DATE_STAMP}.benchmark.txt'
    conda:
      '../envs/label_celltypes.yaml'
    shell: """
      exec &>> {log}
      python3 {scprocess_dir}/scripts/train_xgboost.py aggregate_fulldata \
        {input.pred_fs} \
        --annots_f          {input.annots_f} \
        --subsample_preds_f {input.subsample_preds_f} \
        --output_f          {output.fulldata_f} \
        $( [ "{params.label_map_f}" != "" ] && echo "--label_map_f {params.label_map_f}" )
      """
