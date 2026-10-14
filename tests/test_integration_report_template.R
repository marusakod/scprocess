#!/usr/bin/env Rscript
# Exercise full template substitution, including the score diagnostic input.
source('scripts/render_htmls.R')
args = list(your_name = 'test', affiliation = 'test', short_tag = 'test',
  date_stamp = '2026-04-01', threads = 1L, qc_dt_f = 'qc.csv.gz',
  doublet_real_f = 'scores.csv.gz', integration_f = 'integration.csv.gz',
  int_embedding = 'harmony', int_res_ls_str = '0.2', int_dbl_cl_prop = 0.5,
  demux_type = 'none', batch_var = 'sample_id')
params = do.call(get_sub_ls, c(list(rule = 'integration', proj_dir = tempdir()), args))
stopifnot(identical(params$doublet_real_f, args$doublet_real_f))
rmd = tempfile(fileext = '.Rmd')
make_rmd_from_temp('integration', 'resources/rmd_templates/integration.Rmd.template', params, rmd)
lines = readLines(rmd)
stopifnot(any(grepl('doublet_real_f = "scores.csv.gz"', lines, fixed = TRUE)),
  !any(grepl('${', lines, fixed = TRUE)))
parse(text = paste(lines[grepl('^doublet_real_f', lines)], collapse = '\n'))
message('Full integration template substitution passed, including doublet_real_f.')
