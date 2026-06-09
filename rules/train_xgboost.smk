if TRAIN_XGB_PARAMS is not None:

  rule train_xgboost_train:
    input:
      annots_f    = TRAIN_XGB_PARAMS['annots_f'],
      cluster_csv = f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
      h5ads_yaml  = f'{int_dir}/h5ads_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml',
    output:
      model_f  = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_xgboost_model.json',
      cls_f    = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_allowed_cls.csv',
      genes_f  = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_selected_genes.txt',
      imp_f    = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_gene_importance.csv',
      preds_f  = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_predictions.csv.gz',
      pb_f     = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_pseudobulk.h5ad',
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
      python3 {scprocess_dir}/scripts/train_xgboost.py \
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


  rule render_html_train_xgboost:
    input:
      r_mkr_f       = f"{code_dir}/marker_genes.R",
      r_utils_f     = f"{code_dir}/utils.R",
      preds_f       = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_predictions.csv.gz'
      imp_f         = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_gene_importance.csv',
      pb_f          = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_pseudobulk.h5ad'
    output:
      r_train_f     = f"{code_dir}/train_xgboost.R",
      rmd_f         = f"{rmd_dir}/{SHORT_TAG}_train_xgboost.Rmd",
      html_f        = f"{docs_dir}/{SHORT_TAG}_train_xgboost.html"
    params:
      your_name     = config['project'].get('your_name', ''),
      affiliation   = config['project'].get('affiliation', ''),
      short_tag     = SHORT_TAG,
      ref_tag       = TRAIN_XGB_PARAMS['ref_tag'],
      proj_dir      = PROJ_DIR,
      output_dir    = xgb_dir,
      predictions_f = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_predictions.csv.gz',
      importance_f  = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_gene_importance.csv',
      pseudobulk_f  = f'{xgb_dir}/{TRAIN_XGB_PARAMS["ref_tag"]}_pseudobulk.h5ad',
      integration_f = f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
      has_coarse    = "true" if TRAIN_XGB_PARAMS.get('label_map_f') is not None else "false",
      min_cells     = TRAIN_XGB_PARAMS.get('min_cells_expressed', 10),
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'render_html_train_xgboost', 'memory', attempt),
      runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'render_html_train_xgboost', 'time', attempt)
    log: f'{logs_dir}/train_xgboost/render_html_train_xgboost_{DATE_STAMP}.log'
    benchmark: f'{benchmark_dir}/train_xgboost/render_html_train_xgboost_{DATE_STAMP}.benchmark.txt'
    conda:
      '../envs/rlibs.yaml'
    shell: """
      exec &>> {log}
      # copy R scripts to code directory
      cp scripts/train_xgboost.R {output.r_train_f}

      # define template
      template_f=$(realpath resources/rmd_templates/train_xgboost.Rmd.template)
      rule="train_xgboost"

      # render html
      Rscript --vanilla -e "source('scripts/render_htmls.R'); \\
        render_html(
          rule_name       = '$rule',
          proj_dir        = '{params.proj_dir}',
          temp_f          = '$template_f',
          rmd_f           = '{output.rmd_f}',
          your_name       = '{params.your_name}',
          affiliation     = '{params.affiliation}',
          short_tag       = '{params.short_tag}',
          ref_tag         = '{params.ref_tag}',
          predictions_f   = '{params.predictions_f}',
          importance_f    = '{params.importance_f}',
          pseudobulk_f    = '{params.pseudobulk_f}',
          integration_f   = '{params.integration_f}',
          has_coarse      = '{params.has_coarse}',
          min_cells       = '{params.min_cells}',
          n_cores         = '{threads}',
          scripts_dir     = '{params.output_dir}'
        )"
      """
