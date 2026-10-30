# Persist counts and preparation separately so marker fitting can resume cheaply.
MARKER_EDGER_BP_SCRIPT_F = str(scprocess_dir / 'scripts/marker_genes_edger_bp.R')
MARKER_FILE_TAG = f"{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_{DATE_STAMP}"
MARKER_PB_DIR = f'{mkr_dir}/pb_{MARKER_FILE_TAG}.bpcells'
MARKER_COL_DATA_F = f'{mkr_dir}/pb_prepared_coldata_{MARKER_FILE_TAG}.csv.gz'
MARKER_ROW_DATA_F = f'{mkr_dir}/pb_prepared_rowdata_{MARKER_FILE_TAG}.csv.gz'
MARKER_GENES_F = f'{mkr_dir}/pb_marker_genes_{MARKER_FILE_TAG}.csv.gz'
MARKER_HVGS_F = f'{mkr_dir}/pb_hvgs_{MARKER_FILE_TAG}.csv.gz'
MARKER_PLOT_DATA_F = f'{mkr_dir}/pb_plot_data_{MARKER_FILE_TAG}.csv.gz'


rule make_marker_pseudobulks:
  input:
    h5ads_yaml_f  = f'{int_dir}/h5ads_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml',
    integration_f = f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    marker_script = MARKER_EDGER_BP_SCRIPT_F
  output:
    pb_dir = directory(MARKER_PB_DIR)
  params:
    sel_res     = config['marker_genes']['mkr_sel_res'],
    min_cl_size = config['marker_genes']['mkr_min_cl_size'],
    batch_var   = BATCH_VAR,
    project_id  = FULL_TAG
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'make_marker_pseudobulks', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'make_marker_pseudobulks', 'time', attempt)
  log:
    f"{logs_dir}/marker_genes/make_marker_pseudobulks_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/marker_genes/make_marker_pseudobulks_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/rlibs_bpcells.yaml'
  shell: """
    exec &>> {log}
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    Rscript --vanilla -e "source('{input.marker_script}'); make_pseudobulks_bpcells(
      integration_f = '{input.integration_f}',
      h5ads_yaml_f  = '{input.h5ads_yaml_f}',
      pb_dir        = '{output.pb_dir}',
      sel_res       = '{params.sel_res}',
      min_cl_size   =  {params.min_cl_size},
      batch_var     = '{params.batch_var}',
      project_id    = '{params.project_id}',
      exclude_doublets = TRUE,
      n_cores       =  {threads})"
    """


rule prepare_marker_pseudobulks:
  input:
    pb_dir        = MARKER_PB_DIR,
    marker_script = MARKER_EDGER_BP_SCRIPT_F
  output:
    prepared_coldata_f = MARKER_COL_DATA_F,
    prepared_rowdata_f = MARKER_ROW_DATA_F
  params:
    min_cells        = config['marker_genes']['mkr_min_cells'],
    adjust_project_id = "FALSE"
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'prepare_marker_pseudobulks', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'prepare_marker_pseudobulks', 'time', attempt)
  log:
    f"{logs_dir}/marker_genes/prepare_marker_pseudobulks_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/marker_genes/prepare_marker_pseudobulks_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/rlibs_bpcells.yaml'
  shell: """
    exec &>> {log}
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    Rscript --vanilla -e "source('{input.marker_script}'); prepare_pseudobulks_bpcells(
      pb_dir             = '{input.pb_dir}',
      prepared_coldata_f = '{output.prepared_coldata_f}',
      prepared_rowdata_f = '{output.prepared_rowdata_f}',
      min_cells          =  {params.min_cells},
      adjust_project_id  =  {params.adjust_project_id},
      n_cores            =  {threads})"
    """


rule calc_marker_hvgs:
  input:
    pb_dir             = MARKER_PB_DIR,
    prepared_coldata_f = MARKER_COL_DATA_F,
    marker_script      = MARKER_EDGER_BP_SCRIPT_F
  output:
    pb_hvgs_f = MARKER_HVGS_F
  params:
    gene_info_f = config['mapping_af']['gene_info_f']
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'calc_marker_hvgs', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'calc_marker_hvgs', 'time', attempt)
  log:
    f"{logs_dir}/marker_genes/calc_marker_hvgs_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/marker_genes/calc_marker_hvgs_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/rlibs_bpcells.yaml'
  shell: """
    exec &>> {log}
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    Rscript --vanilla -e "source('{input.marker_script}'); calculate_pseudobulk_hvgs_bpcells(
      pb_dir             = '{input.pb_dir}',
      prepared_coldata_f = '{input.prepared_coldata_f}',
      pb_hvgs_f          = '{output.pb_hvgs_f}',
      biotypes_dt        = load_gene_biotypes_bpcells('{params.gene_info_f}'),
      n_cores            =  {threads})"
    """


rule run_marker_genes:
  input:
    pb_dir             = MARKER_PB_DIR,
    prepared_coldata_f = MARKER_COL_DATA_F,
    prepared_rowdata_f = MARKER_ROW_DATA_F,
    marker_script      = MARKER_EDGER_BP_SCRIPT_F
  output:
    mkrs_f = MARKER_GENES_F
  params:
    gene_info_f       = config['mapping_af']['gene_info_f'],
    adjust_project_id = "FALSE",
    dispersion_use_clusters = str(config['marker_genes']['mkr_dispersion_use_clusters']).upper()
  threads: config['resources']['n_run_marker_genes']
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'run_marker_genes', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'run_marker_genes', 'time', attempt)
  log:
    f"{logs_dir}/marker_genes/run_marker_genes_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/marker_genes/run_marker_genes_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/rlibs_bpcells.yaml'
  shell: """
    exec &>> {log}
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    Rscript --vanilla -e "source('{input.marker_script}'); calculate_marker_genes_bpcells(
      pb_dir             = '{input.pb_dir}',
      prepared_coldata_f = '{input.prepared_coldata_f}',
      prepared_rowdata_f = '{input.prepared_rowdata_f}',
      mkrs_f             = '{output.mkrs_f}',
      biotypes_dt        = load_gene_biotypes_bpcells('{params.gene_info_f}'),
      adjust_project_id  =  {params.adjust_project_id},
      dispersion_use_clusters = {params.dispersion_use_clusters},
      n_cores            =  {threads})"
    """


rule prepare_marker_plot_data:
  input:
    pb_dir             = MARKER_PB_DIR,
    prepared_coldata_f = MARKER_COL_DATA_F,
    mkrs_f             = MARKER_GENES_F,
    pb_hvgs_f          = MARKER_HVGS_F,
    marker_script      = MARKER_EDGER_BP_SCRIPT_F,
    marker_plot_script = str(scprocess_dir / 'scripts/marker_genes.R'),
    utils_script       = str(scprocess_dir / 'scripts/utils.R'),
    custom_marker_fs   = config['marker_genes']['custom_mkr_paths'].split(',') if config['marker_genes']['custom_mkr_paths'] else [],
    ambient_f          = f'{empty_dir}/edger_empty_genes_all_{FULL_TAG}_{DATE_STAMP}.csv.gz' if CAN_CALC_AMBIENT_GENES else []
  output:
    pb_plot_data_f = MARKER_PLOT_DATA_F
  params:
    min_cpm_mkr = config['marker_genes']['mkr_min_cpm_mkr'],
    not_ok_re   = config['marker_genes']['mkr_not_ok_re'],
    custom_mkr_paths = lambda wildcards, input: ','.join(input.custom_marker_fs)
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'prepare_marker_plot_data', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'prepare_marker_plot_data', 'time', attempt)
  log:
    f"{logs_dir}/marker_genes/prepare_marker_plot_data_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/marker_genes/prepare_marker_plot_data_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/rlibs_bpcells.yaml'
  shell: """
    exec &>> {log}
    Rscript --vanilla -e "source('{input.utils_script}'); source('{input.marker_script}'); source('{input.marker_plot_script}'); make_report_logcpms_bpcells(
      pb_dir             = '{input.pb_dir}',
      prepared_coldata_f = '{input.prepared_coldata_f}',
      mkrs_f             = '{input.mkrs_f}',
      pb_hvgs_f          = '{input.pb_hvgs_f}',
      out_f              = '{output.pb_plot_data_f}',
      min_cpm_mkr        = {params.min_cpm_mkr},
      not_ok_re          = '{params.not_ok_re}',
      ambient_f          = '{input.ambient_f}',
      custom_mkr_paths   = strsplit('{params.custom_mkr_paths}', ',', fixed = TRUE)[[1]])"
    """


rule run_fgsea:
  input:
    mkrs_f        = f'{mkr_dir}/pb_marker_genes_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_{DATE_STAMP}.csv.gz'
  output:
    fgsea_go_bp_f = f'{mkr_dir}/fgsea_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_go_bp_{DATE_STAMP}.csv.gz', 
    fgsea_go_cc_f = f'{mkr_dir}/fgsea_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_go_cc_{DATE_STAMP}.csv.gz',
    fgsea_go_mf_f = f'{mkr_dir}/fgsea_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_go_mf_{DATE_STAMP}.csv.gz'
  params:
    genome_ref      = GENOME_REF,
    mkr_gsea_dir    = config['marker_genes']['mkr_gsea_dir'],
    mkr_min_cpm_go  = config['marker_genes']['mkr_min_cpm_go'],
    mkr_max_zero_p  = config['marker_genes']['mkr_max_zero_p'],
    mkr_gsea_var    = config['marker_genes']['mkr_gsea_var'],
    mkr_gsea_cut    = config['marker_genes']['mkr_gsea_cut'], 
    mkr_not_ok_re   = config['marker_genes']['mkr_not_ok_re']
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'run_fgsea', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'run_fgsea', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/marker_genes/run_fgsea_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/marker_genes/run_fgsea_{DATE_STAMP}.log'
  conda: '../envs/rlibs.yaml'
  shell:"""
    exec &>> {log}
    
    Rscript --vanilla -e "source('{scprocess_dir}/scripts/utils.R'); source('{scprocess_dir}/scripts/fgsea.R'); run_fgsea(
      mkrs_f        = '{input.mkrs_f}', 
      fgsea_go_bp_f = '{output.fgsea_go_bp_f}', 
      fgsea_go_cc_f = '{output.fgsea_go_cc_f}', 
      fgsea_go_mf_f = '{output.fgsea_go_mf_f}', 
      genome_ref    = '{params.genome_ref}',
      gsea_dir      = '{params.mkr_gsea_dir}', 
      min_cpm_go    =  {params.mkr_min_cpm_go}, 
      max_zero_p    =  {params.mkr_max_zero_p},
      gsea_var      = '{params.mkr_gsea_var}',
      gsea_cut      =  {params.mkr_gsea_cut},
      not_ok_re     = '{params.mkr_not_ok_re}',
      n_cores       =  {threads})"
    """
