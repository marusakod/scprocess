import os
import pathlib
import yaml

scprocess_dir = pathlib.Path(config.get("scprocess_dir", workflow.basedir)).parent
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
    config_f = config["_config_path"]
  output:
    model_f  = f'{output_dir}/{ref_tag}_xgboost_model.json',
    cls_f    = f'{output_dir}/{ref_tag}_allowed_cls.csv',
    genes_f  = f'{output_dir}/{ref_tag}_selected_genes.txt',
    imp_f    = f'{output_dir}/{ref_tag}_gene_importance.csv',
    preds_f  = f'{output_dir}/{ref_tag}_predictions.csv.gz',
  threads: n_cores
  resources:
    mem_mb  = mem_gb * MB_PER_GB,
    runtime = runtime_mins
  conda:
    '../envs/train_xgboost.yaml'
  shell: """
    python3 {scprocess_dir}/scripts/train_xgboost.py {input.config_f}
    """


rule render_html_train_xgboost:
  input:
    preds_f       = rules.train_xgboost.output.preds_f,
    imp_f         = rules.train_xgboost.output.imp_f,
  output:
    r_utils_f     = f"{code_dir}/utils.R",
    r_lbl_f       = f"{code_dir}/label_celltypes.R",
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
