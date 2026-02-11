# functions for setting up scprocess_data


import argparse
import glob
import gzip
import numpy as np
import os
import pandas as pd
import re
import re
import subprocess
import warnings
import yaml

TENX_NAMES    = ['human_2020', 'human_2024', 'mouse_2020', 'mouse_2024']
TENX_MITOS    = {
  'human_2020': "^MT-", 
  'human_2024': "^MT-",
  'mouse_2020': "^mt-", 
  'mouse_2024': "^mt-"
}
URLS_10X_REFS = {
  'human_2020': "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz", 
  'human_2024': "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz", 
  'mouse_2020': "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz", 
  'mouse_2024': "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
}
URLS_GTF_TXTS = { 
  'human_2020': "https://zenodo.org/records/14247195/files/human_2020_rrna_genes_gtf.txt.gz", 
  'human_2024': "https://zenodo.org/records/14247195/files/human_2024_rrna_genes_gtf.txt.gz", 
  'mouse_2020': "https://zenodo.org/records/14247195/files/mouse_2020_rrna_genes_gtf.txt.gz", 
  'mouse_2024': "https://zenodo.org/records/14247195/files/mouse_2024_rrna_genes_gtf.txt.gz"
}
URLS_ZEN_IDXS = { 
  'human_2020': "https://zenodo.org/records/14247195/files/alevin-idx_human_2020_rrna.tar.gz", 
  'human_2024': "https://zenodo.org/records/14247195/files/alevin-idx_human_2024_rrna.tar.gz", 
  'mouse_2020': "https://zenodo.org/records/14247195/files/alevin-idx_mouse_2020_rrna.tar.gz", 
  'mouse_2024': "https://zenodo.org/records/14247195/files/alevin-idx_mouse_2024_rrna.tar.gz"
}


def get_scprocess_data(scdata_dir):
  print('Downloading data from scprocessData github repo')
  # switch to scprocess data dir
  os.chdir(scdata_dir)

  # download tar archive from scprocessData and extract
  subprocess.run('wget https://github.com/marusakod/scprocessData/releases/download/v0.1.2/scprocess_data_archive.tar.gz',
  shell=True, capture_output=False)
  
  subprocess.run('tar xvf scprocess_data_archive.tar.gz', shell=True, capture_output=False)

  # remove tar archive
  os.remove('scprocess_data_archive.tar.gz')
  
  # check if all necessary directories are there
  dirs_ls = ["cellranger_ref", "gmt_pathways", "marker_genes", "xgboost"]
  for dir in dirs_ls:
    assert os.path.isdir(dir), \
      f"{dir} directory doesn't exist"
  
  print('Done!')

  return


def get_af_index_parameters(config):
  # initialize
  genomes     = config['genomes']
  SETUP_LS    = []

  # get parameters for all specified tenx genomes
  if 'tenx' in genomes:
    # get just tenx list
    tenx_ls     = genomes['tenx']
    for spec_tenx in tenx_ls:
      SETUP_LS.append(_get_index_parameters_tenx(spec_tenx, TENX_NAMES))

  # get parameters for all specified custom genomes
  if 'custom' in genomes:
    # get just custom list
    custom_ls     = genomes['custom']
    for spec_custom in custom_ls:
      SETUP_LS.append(_get_index_parameters_custom(spec_custom, TENX_NAMES))

  # check no duplicate names
  setup_names = [ s['name'] for s in SETUP_LS ]
  if not len(setup_names) == len(set(setup_names)):
    raise KeyError("Duplicated genome names are not allowed!")

  # check if reference for tutorial needs to be added
  if 'mouse_2024' not in setup_names:
    spec_tmp    = {'name': "mouse_2024"}
    setup_names = [ setup_names, 'mouse_2024' ]
    SETUP_LS.append( _get_index_parameters_tenx(spec_tmp, TENX_NAMES) )

  # make dictionary instead
  SETUP_LS    = dict(zip(setup_names, SETUP_LS))

  return SETUP_LS


def _get_index_parameters_tenx(spec_tenx, TENX_NAMES):
  # define defaults
  spec_tenx["gtf"]          = None
  spec_tenx["fasta"]        = None
  spec_tenx["mito_str"]     = TENX_MITOS[ spec_tenx['name'] ]
  spec_tenx["is_tenx"]      = True
  if not "decoys" in spec_tenx:
    spec_tenx["decoys"]       = True
  if not "rrnas" in spec_tenx:
    spec_tenx["rrnas"]        = True
  spec_tenx['is_prebuilt']  = spec_tenx['decoys'] and spec_tenx['rrnas']

  return spec_tenx


def _get_index_parameters_custom(spec_custom, TENX_NAMES):
  # check if all required entries are specified
  if not os.path.isfile(spec_custom['gtf']):
    raise FileNotFoundError(f"file {gtf} specified in configfile doesn't exist")
  if not os.path.isfile(spec_custom['fasta']):
    raise FileNotFoundError(f"file {fa} specified in configfile doesn't exist")

  # fill in defaults if we need to
  if not "decoys" in spec_custom:
    spec_custom["decoys"]     = True
  if not "rrnas" in spec_custom:
    spec_custom["rrnas"]      = True

  return spec_custom


def _safe_boolean(val):
  if type(val) is bool:
    res = val
  elif val in ["True", "true"]:
    res = True
  elif val in ["False", "false"]:
    res = False
  else:
    raise ValueError('{val} is not a boolean')

  return res


def _check_valid_index(idx_path):
  # check if main directory exists
  if not os.path.isdir(idx_path):
    raise ValueError(f"The provided path '{idx_path}' is not a valid directory.")

  # define expected directories and files
  req_dirs = {
    "idx_dir":  os.path.join(idx_path, 'index'),
    "ref_dir":  os.path.join(idx_path, 'ref')
  }
  req_fs = {
    "index_info.json":          os.path.join(idx_path, 'index_info.json'),
    "simpleaf_index_log.json":  os.path.join(idx_path, 'simpleaf_index_log.json'),
    "gene_id_to_name.tsv":      os.path.join(req_dirs["ref_dir"], 'gene_id_to_name.tsv'),
    "roers_make-ref.json":      os.path.join(req_dirs["ref_dir"], 'roers_make-ref.json'),
    "roers_ref.fa":             os.path.join(req_dirs["ref_dir"], 'roers_ref.fa'),
    "t2g_3col.tsv":             os.path.join(req_dirs["ref_dir"], 't2g_3col.tsv')
  }

  # Check required directories
  for dir_name, dir_path in req_dirs.items():
    if not os.path.isdir(dir_path):
      raise ValueError(f"The required directory '{dir_name}' is missing at '{idx_path}'.")

  # Check required files
  missing_fs = [file_name for file_name, file_path in req_fs.items() if not os.path.isfile(file_path)]
  if missing_fs:
    raise ValueError(f"The following required files are missing: {', '.join(missing_fs)}")

  return True


# function that makes simpleaf index
def set_up_af_index(scdata_dir, txome_name, fasta_f, gtf_f, index_dir, mito_str, is_prebuilt, is_tenx, has_decoy, has_rrna, n_cores):
  if is_prebuilt:
    print(f"'is_prebuilt' is True, value is {is_prebuilt}")
  else:
    print(f"'is_prebuilt' is False, value is {is_prebuilt}")

  # create output directories
  ref_dir   = os.path.join(scdata_dir, 'reference_genomes', txome_name)
  idx_dir   = os.path.join(scdata_dir, 'alevin_fry_home', txome_name)
  os.makedirs(ref_dir,  exist_ok=True)
  os.makedirs(idx_dir,  exist_ok=True)
  gtf_txt_f = os.path.join(ref_dir, f'{txome_name}_genes_gtf.txt.gz')

  # if prebuilt, just download index
  if is_prebuilt:
    _download_prebuilt_index(txome_name, idx_dir)
    _download_gtf_txt_file(txome_name, ref_dir, gtf_txt_f)
  else:
    if is_tenx:
      if has_rrna == True:
        fasta_f, gtf_f  = _make_10x_fasta_and_gtf_w_rrna(ref_dir, txome_name)
      else:
        fasta_f, gtf_f  = _download_predefined_fasta_and_gtf(ref_dir, txome_name)
    else:
      # copy custom files
      fasta_f, gtf_f  = _copy_custom_fasta_and_gtf(ref_dir, txome_name, fasta_f, gtf_f)

    # build index
    _build_index_w_simpleaf(txome_name, idx_dir, has_decoy, fasta_f, gtf_f, n_cores = n_cores)

    # save gtf_txt file
    _make_gtf_txt_file(gtf_f, gtf_txt_f)


  # create one dataframe from all dictionaries
  yaml_f    = os.path.join(idx_dir, f"{txome_name}_index_params.yaml")
  _make_index_params_yaml(yaml_f, txome_name, fasta_f, index_dir, gtf_f, gtf_txt_f, mito_str, has_decoy, has_rrna, is_prebuilt, is_tenx)

  print(f'completed making index for {txome_name} genome in {scdata_dir}.')

  return 


def _download_prebuilt_index(ref_txome, idx_dir):
  # get index url
  print('Downloading alevin index for ' + ref_txome)
  idx_url   = URLS_ZEN_IDXS[ref_txome]

  # make output directory
  os.chdir(idx_dir)

  # download index
  subprocess.run(f"wget {idx_url}", shell=True)
  idx_name  = f'alevin-idx_{ref_txome}_rrna.tar.gz'

  # untar
  subprocess.run(f'tar --strip-components=1 -xvf {idx_name}', shell=True, capture_output=False)

  # remove tar archive
  os.remove(idx_name)

  return


def _download_gtf_txt_file(ref_txome, ref_dir, gtf_txt_f):
  # get index url
  print('Downloading GTF txt file for ' + ref_txome)
  txt_url   = URLS_GTF_TXTS[ref_txome]

  # make output directory
  os.chdir(ref_dir)
  subprocess.run(f"wget -O {gtf_txt_f} {txt_url}", shell=True)

  return


def _make_10x_fasta_and_gtf_w_rrna(ref_dir, ref_txome):
  print(f"Creating {ref_txome} 10x genome with rRNAs")
  # get bash script to download gtf and fasta
  bash_f     = f"./scripts/build_10x_style_genomes/build_10x_style_{ref_txome}_genome_w_rRNAs.sh"
  bash_f     = os.path.realpath(bash_f)

  # create a new directory for the specified genome and switch to it
  gnome_dir  = os.path.join(ref_dir, ref_txome)
  os.makedirs(gnome_dir, exist_ok=True)
  os.chdir(gnome_dir)

  # run bash script that downloads gtf and 
  subprocess.run(bash_f, shell=True)

  # define files
  fasta_f   = os.path.join(gnome_dir, 'genome.fa')
  gtf_f     = os.path.join(gnome_dir, 'genes.gtf')
  print('Done!')

  return fasta_f, gtf_f


def _download_predefined_fasta_and_gtf(ref_dir, ref_txome):
  print(f'Downloading {ref_txome} genome from 10x')
  link          = URLS_10X_REFS[ref_txome]
  # get tarball 
  tball         = link.rsplit('/', 1)[-1]
  
  # create a new directory for the specified genome and switch to it
  gnome_dir     = os.path.join(ref_dir, ref_txome)
  os.makedirs(gnome_dir, exist_ok = True)
  os.chdir(gnome_dir)
        
  # download reference into that dir
  subprocess.run(f'wget {link}', shell=True)
        
  # list tar contents
  tar_out       = subprocess.run(f'tar -tzf {tball}', shell=True, capture_output=True, text=True)
  tball_all_fs  = tar_out.stdout.splitlines()
        
  # look for genome.fa and genes.gtf (or genes.gtf.gz) files and extract the
  patt          = re.compile(r'(genome\.fa|genes\.gtf(\.gz)?)$') 
  tball_filt_fs = [file for file in tball_all_fs if patt.search(file)]

  # extract everything
  for t in tball_filt_fs:
    ext_spell     = ['tar', '--extract', '--file=' + tball, '--strip-components=2', t]
    subprocess.run(ext_spell, check=True)

  # unzip
  gnome_fs = os.listdir()
  for f in gnome_fs:
    if f.endswith('.gtf.gz'):
      subprocess.run(f'gunzip {f}', shell=True, capture_output=False)

  # delete the tarball
  os.remove(tball)
  print('Done!')

  fasta_f   = os.path.join(gnome_dir, 'genome.fa')
  gtf_f     = os.path.join(gnome_dir, 'genes.gtf')

  return fasta_f, gtf_f


def _copy_custom_fasta_and_gtf(ref_dir, ref_txome, fasta_f, gtf_f):
  # get index url
  print('Copying fasta and gtf for ' + ref_txome)

  # create a new directory for the specified genome and switch to it
  os.makedirs(ref_dir, exist_ok = True)
  os.chdir(ref_dir)
  subprocess.run(f"rsync {fasta_f} genome.fa", shell=True)
  subprocess.run(f"rsync {gtf_f} genes.gtf", shell=True)
  print('Done.')

  fasta_f   = os.path.join(ref_dir, 'genome.fa')
  gtf_f     = os.path.join(ref_dir, 'genes.gtf')

  return fasta_f, gtf_f


def _build_index_w_simpleaf(ref_txome, idx_dir, has_decoy, fasta_f, gtf_f, n_cores = 8):
  print(f'Creating alevin index for {ref_txome} { "with decoys " if has_decoy else "" }in {idx_dir}')
   
  # define whether or not to include --decoy-paths flag
  decoy_flag  = f"--decoy-paths {fasta_f} " if has_decoy else ""

  # code for making index
  bash_script = f"""
  #!/bin/bash
  ulimit -n 2048

  # simpleaf configuration
  export ALEVIN_FRY_HOME="/tmp/alevin_fry_home"
  mkdir -p ${{ALEVIN_FRY_HOME}}
  if [ ! -d "${{ALEVIN_FRY_HOME}}" ]; then
    raise error "couldn't create ALEVIN_FRY_HOME directory in /tmp"
  fi
  simpleaf set-paths

  # change working directory to tmp directory
  cd ${{ALEVIN_FRY_HOME}}

  # set up this build
  TMP_IDX_DIR="${{ALEVIN_FRY_HOME}}/{genome_name}"

  # simpleaf index
  simpleaf index \
    --output ${{TMP_IDX_DIR}} \
    --fasta {fasta_f} \
    --gtf {gtf_f} \
    {decoy_flag} --overwrite --threads {n_cores} \
    --use-piscem

  # copy results to nice place, tidy up
  rsync -avP ${{TMP_IDX_DIR}}/ {idx_dir}
  rm -rf ${{ALEVIN_FRY_HOME}}
  """

  # run bash script
  print(bash_script)
  subprocess.run(bash_script, shell=True, executable='/bin/bash', check=True)
  
  return


def _make_gtf_txt_file(gtf_f, gtf_txt_f):
  print(f'Creating txt files from gtf')
  # load gtf
  import pyranges as pr
  gtf         = pr.read_gtf(gtf_f, as_df = True)

  # restrict to just genes, remove some columns
  gene_annots = gtf[ gtf['Feature'] == 'gene' ]
  if "gene_biotype" in gene_annots.columns:
    gene_annots.rename(columns={
      'gene_biotype': 'gene_type'
    }, inplace=True)

  gene_annots = gene_annots[['gene_id', 'gene_name', 'gene_type', 'Chromosome', 'Start', 'End', 'Strand']]
  
  # Rename columns
  gene_annots.rename(columns={
    'gene_id': 'ensembl_id',
    'gene_name': 'symbol',
    'Chromosome': 'chromosome',
    'Start': 'start',
    'End': 'end',
    'Strand': 'strand'
  }, inplace=True)

  # check for any NAs in symbol
  gene_annots["symbol"] = gene_annots["symbol"].fillna(gene_annots["ensembl_id"])

  # calculate width and add gene_id column
  gene_annots['width']    = gene_annots['end'] - gene_annots['start'] + 1
  gene_annots['gene_id']  = gene_annots['symbol'] + '_' + gene_annots['ensembl_id']
  
  # do some sorting of chromosomes 
  chr_freq = gene_annots['chromosome'].value_counts()
  gene_annots['chr_freq'] = gene_annots['chromosome'].map(chr_freq)
  gene_annots_sorted = gene_annots.sort_values(by=['chr_freq', 'start', 'end'], ascending=[False, True, True])
  gene_annots_sorted = gene_annots_sorted.drop(columns=['chr_freq'])
  
  # make gene_id the first column
  gene_id_col = gene_annots_sorted.pop('gene_id')
  gene_annots_sorted.insert(0, 'gene_id', gene_id_col)
  
  # make sure there are no duplicated gene ids
  assert gene_annots_sorted['gene_id'].duplicated().sum() == 0

  # save
  with gzip.open(gtf_txt_f, 'wt') as f:
    gene_annots_sorted.to_csv(f, sep='\t', index=False)
      
  return


def _make_index_params_yaml(yaml_f, ref_txome, fasta_f, index_dir, gtf_f, gtf_txt_f, mito_str, 
  has_decoy, has_rrna, is_prebuilt, is_tenx):
  param_ls  = {
    "ref_txome":    ref_txome, 
    "fasta_f":      fasta_f, 
    "index_dir":    index_dir, 
    "gtf_f":        gtf_f, 
    "gtf_txt_f":    gtf_txt_f, 
    "mito_str":     mito_str, 
    "has_decoy":    has_decoy, 
    "has_rrna":     has_rrna, 
    "is_prebuilt":  is_prebuilt, 
    "is_tenx":      is_tenx
  }
  with open(yaml_f, 'w') as f:
    yaml.dump(param_ls, f)

  return


def save_index_params_csv(csv_f, yaml_fs):
  # load all yamls
  assert all([ os.path.isfile(f) for f in yaml_fs]), "not all yaml files exist"

  # make df
  params_df = pd.DataFrame( columns = ["ref_txome", "mito_str", "gtf_txt_f"] )
  for i, yaml_f in enumerate(yaml_fs):
    with open(yaml_f) as f:
      yaml_ls = yaml.load(f, Loader=yaml.FullLoader)
    params_df.loc[i] = [yaml_ls['ref_txome'], yaml_ls['mito_str'], yaml_ls['gtf_txt_f']]

  # save to csv
  params_df.to_csv(csv_f, index = False)


if __name__ == "__main__":
  # set up parser
  parser      = argparse.ArgumentParser()
  subparsers  = parser.add_subparsers(dest="function_name", help="Name of the function to run")

  # parser for get_scprocess_data
  getdata     = subparsers.add_parser('get_scprocess_data')
  getdata.add_argument('scdata_dir', type=str)

  # parsers for set_up_af_index
  get_af      = subparsers.add_parser('set_up_af_index')
  get_af.add_argument('scdata_dir', type = str)
  get_af.add_argument('genome', type = str, help = 'genome name')
  get_af.add_argument('fasta_f', type = str, help = 'path to fasta file')
  get_af.add_argument('gtf_f', type = str, help = 'path to gtf file')
  get_af.add_argument('index_dir', type =str, help = 'path to prebuilt alevin index' )
  get_af.add_argument('mito_str', type = str, help = 'mitochondrial gene identifier')
  get_af.add_argument('is_prebuilt', type = str, 
    help = 'bool value defining whether or not the user has specified a prebuilt index')
  get_af.add_argument('is_tenx', type = str, 
    help = 'bool value defining whether or not the user has specified an ref genome provided by 10x')
  get_af.add_argument('has_decoy', type = str, 
    help = 'bool value defining whether or not to use decoys when building indices with simpleaf')
  get_af.add_argument('has_rrna', type = str, 
    help = 'bool value defining whether or not to include ribosomal rrnas for simpleaf index')
  get_af.add_argument('cores', type = int)

  # parsers for save_index_params_csv
  save_csv    = subparsers.add_parser('save_index_params_csv')
  save_csv.add_argument('csv_f', type = str)
  save_csv.add_argument('yaml_fs', type = str, nargs = "+", help = 'list of yaml files')

  # decide which function
  args = parser.parse_args()
  if args.function_name == 'get_scprocess_data':
    get_scprocess_data(args.scdata_dir)
  elif args.function_name == 'set_up_af_index':
    set_up_af_index(args.scdata_dir, args.genome, args.fasta_f, args.gtf_f, args.index_dir, args.mito_str, 
      _safe_boolean(args.is_prebuilt), _safe_boolean(args.is_tenx), 
      _safe_boolean(args.has_decoy), _safe_boolean(args.has_rrna), args.cores)
  elif args.function_name == 'save_index_params_csv':
    save_index_params_csv(args.csv_f, args.yaml_fs)
  else:
    parser.print_help()

