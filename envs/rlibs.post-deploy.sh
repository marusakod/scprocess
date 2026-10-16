#!/usr/bin/env bash
# Snakemake activates the newly created environment before running this hook.
set -euo pipefail
: "${CONDA_PREFIX:?Run this hook inside the new Snakemake Conda environment}"
scdbl_commit=e66180a1d84733bbcf63ed8590ae75a4a5329e34
source_dir=$(mktemp -d "${TMPDIR:-/tmp}/scprocess-scdblfinder.XXXXXX")
trap 'rm -rf "$source_dir"' EXIT
# Prevent shared/user R libraries from replacing the environment's packages.
export R_LIBS="" R_LIBS_USER="" R_LIBS_SITE="" R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset R_HOME
git clone --quiet --filter=blob:none https://github.com/plger/scDblFinder.git "$source_dir"
git -C "$source_dir" checkout --quiet --detach "$scdbl_commit"
test "$(git -C "$source_dir" rev-parse HEAD)" = "$scdbl_commit"
printf 'RemoteType: github\nRemoteUsername: plger\nRemoteRepo: scDblFinder\nRemoteSha: %s\n' "$scdbl_commit" >> "$source_dir/DESCRIPTION"
"$CONDA_PREFIX/bin/R" CMD INSTALL --no-multiarch --library="$CONDA_PREFIX/lib/R/library" "$source_dir"
export SCPROCESS_SCDBL_COMMIT="$scdbl_commit"
"$CONDA_PREFIX/bin/Rscript" --vanilla - <<'RSCRIPT'
suppressPackageStartupMessages(library(scDblFinder))
expected = c(scDblFinder = '1.27.12', xgboost = '3.4.2.1')
observed = vapply(names(expected), function(p) as.character(packageVersion(p)), character(1))
if (!identical(observed, expected)) stop('Doublet package version mismatch: ', paste(observed, collapse = ', '))
commit = Sys.getenv('SCPROCESS_SCDBL_COMMIT')
if (!identical(packageDescription('scDblFinder')$RemoteSha, commit)) stop('scDblFinder GitHub commit mismatch')
# Test both fixed-round and CV training directly. scDblFinder can otherwise
# swallow a training error and return its preliminary score without warning.
set.seed(219182L)
x = matrix(rnorm(1000), nrow = 200)
y = factor(ifelse(x[,1] > 0, 'doublet', 'real'), levels = c('real','doublet'))
train = get('.xgbtrain', asNamespace('scDblFinder'))
for (rounds in list(10L, NULL)) {
  set.seed(219182L)
  fit = train(x, y, nrounds = rounds, nthreads = 1)
  scores = predict(fit, x)
  stopifnot(length(scores) == nrow(x), all(is.finite(scores)), all(scores >= 0 & scores <= 1))
}
manifest_dir = file.path(Sys.getenv('CONDA_PREFIX'), 'share/scprocess')
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(package = names(observed), version = unname(observed),
  source = c(paste0('https://github.com/plger/scDblFinder/commit/', commit), 'conda-forge/r-xgboost=3.4.2')),
  file.path(manifest_dir, 'doublet-package-versions.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
message('Verified scDblFinder ', observed[['scDblFinder']], ' at ', commit,
  ' with XGBoost ', observed[['xgboost']], '; both classifier smoke tests passed.')
RSCRIPT
