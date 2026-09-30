# Tiny design check; no expression analysis or package installation.
expressions <- parse('scripts/marker_genes_edger_bp.R')
definition <- Filter(function(x) is.call(x) && identical(x[[1]], as.name('<-')) &&
  identical(x[[2]], as.name('.join_project_effect')), as.list(expressions))
stopifnot(length(definition) == 1L)
eval(definition[[1]])
stopifnot(identical(.join_project_effect(TRUE, c('a', 'a')), FALSE))
stopifnot(identical(.join_project_effect(TRUE, c('a', 'b')), TRUE))
stopifnot(identical(.join_project_effect(FALSE, c('a', 'b')), FALSE))
group <- factor(c('cl1', 'cl1', 'cl2', 'cl2'))
project_id <- factor(rep('a', 4))
design <- if (.join_project_effect(TRUE, project_id)) model.matrix(~project_id + group) else model.matrix(~group)
stopifnot(qr(design)$rank == ncol(design), ncol(design) == 2L)
cat('Four joined-zoom project-effect checks passed.\n')
