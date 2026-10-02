# Render the template against small existing test outputs in an isolated project.
# Usage: Rscript --vanilla tests/test_cell_cycle_report.R SOURCE_PROJECT TEST_PROJECT
args = commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
repo = normalizePath(".")
source_project = normalizePath(args[1])
test_project = args[2]
stopifnot(!dir.exists(test_project))
source(file.path(repo, "scripts/cell_cycle.R"))
workflowr::wflow_start(test_project, change_wd = FALSE,
  user.name = "scprocess test", user.email = "test@example.invalid")
test_project = normalizePath(test_project)
site_f = file.path(test_project, "analysis/_site.yml")
site = yaml::read_yaml(site_f)
site$output_dir = "../public"
yaml::write_yaml(site, site_f)
file.copy(file.path(repo, "scripts", c("utils.R", "cell_cycle.R")),
          file.path(test_project, "code"))
cycle_dir = file.path(source_project, "output/test_cell_cycle")
tag = "test_project_2026-01-01"
# Build a small report fixture from saved coordinates; source outputs stay untouched.
blob_f = file.path(test_project, 'data/blob.csv')
blob = estimate_cell_cycle_blob(
  file.path(cycle_dir, paste0('tricycle_origin_diagnostics_', tag, '.csv.gz')), blob_f)
diagnostics = fread(file.path(cycle_dir,
  paste0('tricycle_origin_diagnostics_', tag, '.csv.gz')))
# Constant inclusion probabilities give an ordinary mean; use a small target
# to catch percentage captions that incorrectly round 0.005 to zero percent.
# This tests report rendering; the Python tests exercise the KDE estimator.
diagnostics[, `:=`(kde_density=1, inclusion_probability=0.005, diagnostic_retained=FALSE)]
diagnostics[seq_len(round(.N * 0.005)), diagnostic_retained := TRUE]
centre = colMeans(as.matrix(diagnostics[, .(tricycle_pc1, tricycle_pc2)]))
origin = data.table(origin_method='kde_density_equalized', target_prop=0.005, seed=20230308L,
  expected_retained_cells=nrow(diagnostics)*0.005,
  tricycle_pc1_origin=centre[1], tricycle_pc2_origin=centre[2])
sensitivity = data.table(target_prop=c(0.005, 1),
  expected_retained_cells=nrow(diagnostics)*c(0.005, 1),
  tricycle_pc1_origin=centre[1], tricycle_pc2_origin=centre[2],
  selected=c(TRUE, FALSE))
diagnostics_f = file.path(test_project, 'data/diagnostics.csv.gz')
fwrite(diagnostics, diagnostics_f)
scores = fread(file.path(cycle_dir, paste0('tricycle_scores_', tag, '.csv.gz')))
scores[, tricycle_theta := atan2(tricycle_pc2-centre[2], tricycle_pc1-centre[1]) %% (2*pi)]
origin_f = file.path(test_project, 'data/origin.csv')
sensitivity_f = file.path(test_project, 'data/sensitivity.csv')
scores_f = file.path(test_project, 'data/scores.csv.gz')
fwrite(origin, origin_f)
fwrite(sensitivity, sensitivity_f)
fwrite(scores, scores_f)
source(file.path(repo, "scripts/render_htmls.R"))
render_html(
  rule_name = "cell_cycle", proj_dir = test_project,
  temp_f = file.path(repo, "resources/rmd_templates/cell_cycle.Rmd.template"),
  rmd_f = file.path(test_project, "analysis/test_cell_cycle.Rmd"),
  your_name = "scprocess test", affiliation = "Validation", short_tag = "test",
  date_stamp = "2026-01-01", threads = 1,
  scores_f = scores_f,
  marker_expression_f = file.path(cycle_dir, paste0("tricycle_marker_expression_", tag, ".csv.gz")),
  diagnostics_f = diagnostics_f,
  blob_f = blob_f,
  origin_f = origin_f, sensitivity_f = sensitivity_f,
  coldata_f = file.path(source_project, "output/test_qc", paste0("coldata_dt_all_cells_", tag, ".csv.gz")))
html = readLines(file.path(test_project, "public/test_cell_cycle.html"), warn = FALSE)
stopifnot(any(grepl("</html>", html, fixed = TRUE)))
document = xml2::read_html(file.path(test_project, "public/test_cell_cycle.html"))
images = xml2::xml_attr(xml2::xml_find_all(document, ".//img"), "src")
images = images[!grepl("^(https?:|data:)", images)]
stopifnot(length(images) > 0, all(file.exists(file.path(test_project, "public", images))))
stopifnot(file.exists(file.path(test_project,
  "public/figure/test_cell_cycle.Rmd/pooled_projection-1.png")))
cat("Cell-cycle origin report render passed:", test_project, "\n")
