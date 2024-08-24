# rule to aggregate large sce objects as metacells

# load modules
import yaml
import pandas as pd
import os
import re
import glob
from snakemake.utils import validate, min_version

# do labelling
rule meta_save_metacells:
  input:
    sce_clean_f = (int_dir + '/sce_clean_' + FULL_TAG + '_' + DATE_STAMP + '.rds') if ('{celltype}' == 'all') else \
      (int_dir + '/sce_clean_' + FULL_TAG + '_' + DATE_STAMP + '.rds')
  output:
    metacell_f  = meta_dir + '/metacells_sce_' + FULL_TAG + '_{celltype}_{max_cells}_' + DATE_STAMP + '.rds',
    mc_map_f    = meta_dir + '/metacells_map_' + FULL_TAG + '_{celltype}_{max_cells}_' + DATE_STAMP + '.txt.gz'
  params:
    celltype    = '{celltype}',
    max_cells   = '{max_cells}'
  threads: 4
  conda:
    '../envs/rlibs.yml'
  resources:
    mem_mb      = 8192
  shell:
    """
    # save sce object
    Rscript -e "\
    source('scripts/utils.R');\
    source('scripts/integration.R');\
    source('scripts/metacells.R');\
    apply_supercell('{input.sce_clean_f}', '{output.metacell_f}', '{output.mc_map_f}',\
      {params.max_cells}, n_cores = {threads})"
    """

