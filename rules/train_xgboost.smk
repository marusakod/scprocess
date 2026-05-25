import os
import pathlib
import yaml

scprocess_dir = pathlib.Path(config["scprocess_dir"]) if "scprocess_dir" in config else pathlib.Path(workflow.basedir).parent
ref_tag       = config["ref_tag"]
output_dir    = config["output_dir"]
n_cores       = config.get("n_cores", 16)

# scprocess project paths (for the html report)
scp_proj      = config.get("_scprocess_project", {})
proj_dir      = scp_proj.get("proj_dir", "")
short_tag     = scp_proj.get("short_tag", "")
rmd_dir       = f"{proj_dir}/analysis"
docs_dir      = f"{proj_dir}/public"
code_dir      = f"{proj_dir}/code"

# resource defaults (from schema or user override)
resources_cfg = config.get("resources", {})
MB_PER_GB     = 1024
mem_gb        = resources_cfg.get("gb_train_xgboost", 64)
runtime_mins  = resources_cfg.get("mins_train_xgboost", 240)

# whether coarse labels are available
has_coarse    = "true" if config.get("label_map_f") is not None else "false"


rule train_xgboost:
  input:
    annots_f    = config["annots_f"],
    cluster_csv = config["cluster_csv"],
    h5ads_yaml  = config["h5ads_yaml"],
  output:
    model_f  = f'{output_dir}/{ref_tag}_xgboost_model.json',
    cls_f    = f'{output_dir}/{ref_tag}_allowed_cls.csv',
    genes_f  = f'{output_dir}/{ref_tag}_selected_genes.txt',
    imp_f    = f'{output_dir}/{ref_tag}_gene_importance.csv',
    preds_f  = f'{output_dir}/{ref_tag}_predictions.csv.gz',
  params:
    ref_tag              = ref_tag,
    output_dir           = output_dir,
    batch_var            = config["batch_var"],
    int_res_ls           = config["int_res_ls"],
    label_map_f          = config.get("label_map_f", ""),
    refine_labels        = config["refine_labels"],
    purity_threshold     = config["purity_threshold"],
    n_cells_per_type     = config["n_cells_per_type"],
    min_cells_per_type   = config["min_cells_per_type"],
    min_cells_expressed  = config["min_cells_expressed"],
    seed                 = config["seed"],
    use_gpu              = config["use_gpu"],
    pass1_subsample      = config["pass1_subsample"],
    pass1_colsample_bytree = config["pass1_colsample_bytree"],
    pass1_learning_rate  = config["pass1_learning_rate"],
    pass1_nrounds        = config["pass1_nrounds"],
    pass1_early_stopping = config["pass1_early_stopping"],
    pass2_colsample_bytree = config["pass2_colsample_bytree"],
    pass2_learning_rate  = config["pass2_learning_rate"],
    pass2_nrounds        = config["pass2_nrounds"],
    pass2_early_stopping = config["pass2_early_stopping"],
    gain_threshold       = config["gain_threshold"],
    min_genes            = config["min_genes"],
    max_genes            = config["max_genes"],
  threads: n_cores
  resources:
    mem_mb  = mem_gb * MB_PER_GB,
    runtime = 180
  conda:
    '../envs/train_xgboost.yaml'
  shell: """
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
    preds_f       = rules.train_xgboost.output.preds_f,
    imp_f         = rules.train_xgboost.output.imp_f,
  output:
    r_utils_f     = f"{code_dir}/utils_tmp.R",
    r_lbl_f       = f"{code_dir}/label_celltypes_tmp.R",
    rmd_f         = f"{rmd_dir}/{short_tag}_train_xgboost.Rmd",
    html_f        = f"{docs_dir}/{short_tag}_train_xgboost.html"
  params:
    your_name     = scp_proj.get("your_name", ""),
    affiliation   = scp_proj.get("affiliation", ""),
    short_tag     = short_tag,
    ref_tag       = ref_tag,
    proj_dir      = proj_dir,
    predictions_f = f'{output_dir}/{ref_tag}_predictions.csv.gz',
    importance_f  = f'{output_dir}/{ref_tag}_gene_importance.csv',
    integration_f = scp_proj.get("cluster_csv", ""),
    has_coarse    = has_coarse
  threads: 1
  resources:
    mem_mb  = 16 * MB_PER_GB,
    runtime = 30
  conda:
    '../envs/rlibs.yaml'
  shell: """
    # copy R code needed by the template
    cp {scprocess_dir}/scripts/utils.R {output.r_utils_f}
    cp {scprocess_dir}/scripts/label_celltypes.R {output.r_lbl_f}

    # define template
    template_f=$(realpath {scprocess_dir}/resources/rmd_templates/train_xgboost.Rmd.template)
    rule="train_xgboost"

    # render html
    Rscript --vanilla -e "source('{scprocess_dir}/scripts/render_htmls.R'); \\
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
        integration_f   = '{params.integration_f}',
        has_coarse      = '{params.has_coarse}'
      )"
    """


rule all:
  input:
    rules.train_xgboost.output.model_f,
    rules.render_html_train_xgboost.output.html_f
  default_target: True
