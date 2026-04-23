# functions for setting up scprocess_data


import argparse
import gzip
import os
import re
import subprocess
import yaml
import tarfile
import shutil
import gzip
import polars as pl

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

TENX_PROBE_SET_NAMES = ['human_v1', 'mouse_v1', 'human_v2', 'mouse_v2']

URLS_PROBE_GTF_TXTS = {
  'human_v1': "https://zenodo.org/records/14247195/files/human_2024_rrna_genes_gtf.txt.gz",
  'human_v2': "https://zenodo.org/records/14247195/files/human_2024_rrna_genes_gtf.txt.gz",
  'mouse_v1': "https://zenodo.org/records/14247195/files/mouse_2024_rrna_genes_gtf.txt.gz",
  'mouse_v2': "https://zenodo.org/records/14247195/files/mouse_2024_rrna_genes_gtf.txt.gz"
}

TENX_PROBE_SET_MITOS = {
  'human_v1': "^MT-",
  'human_v2': "^MT-",
  'mouse_v1': "^mt-",
  'mouse_v2': "^mt-"
}
URLS_10X_PROBE_SETS = {
  'human_v1': "https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv",
  'mouse_v1': "https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Mouse_Transcriptome_Probe_Set_v1.1.1_GRCm39-2024-A.csv",
  'human_v2': "https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv",
  'mouse_v2': "https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Mouse_Transcriptome_Probe_Set_v2.0.0_GRCm39-2024-A.csv"
}


def get_scprocess_data(scdata_dir, ranger_url, whitelists_lu_f, ranger_version_f):
  print('Downloading data from scprocessData github repo')

  if ranger_url == "":
    raise ValueError('Something went wrong. Please provide a valid cellranger link')

  # switch to scprocess data dir
  os.chdir(scdata_dir)

  # download tar archive from scprocessData and extract
  subprocess.run('wget https://github.com/marusakod/scprocessData/releases/download/v0.1.0/scprocess_data_archive.tar.gz',
  shell=True, capture_output=False)
  
  subprocess.run('tar xvf scprocess_data_archive.tar.gz', shell=True, capture_output=False)

  # remove tar archive
  os.remove('scprocess_data_archive.tar.gz')
  
  # check if all necessary directories are there
  dirs_ls = [ "gmt_pathways", "marker_genes", "xgboost"]
  for dir in dirs_ls:
    assert os.path.isdir(dir), \
      f"{dir} directory doesn't exist"
  
  # get version of cellranger
  ranger_version = re.search(r"cellranger-(\d+\.\d+\.\d+)", ranger_url).group(1)
  
  # download cellranger and extract whitelists
  ranger_data_dir = os.path.join(scdata_dir, 'cellranger_ref')
  get_cellranger_whitelists(ranger_data_dir, whitelists_lu_f, ranger_url, ranger_version)

  # save file with cellranger version
  with open(ranger_version_f, "w") as f:
    f.write(str(ranger_version))

  print('Done!')

  return


def get_cellranger_whitelists(output_dir, whitelists_lu_f, ranger_url, ranger_version):

  os.makedirs(output_dir, exist_ok=True)

  tmp_tar    = f"cellranger-{ranger_version}.tar.gz"
  tar_path   = os.path.join(output_dir, tmp_tar)

  # download cellranger
  print(f"Downloading CellRanger")
  subprocess.run(["wget", "-O", tar_path, ranger_url])
  
   # whitelist names in cellranger and new names for scprocess
  cr_wls_dict = {
    "3v2_5v1_5v2": "737K-august-2016.txt",          
    "3v3": "3M-february-2018_TRU.txt.gz",          
    "3v4": "3M-3pgex-may-2023_TRU.txt.gz",         
    "5v3": "3M-5pgex-jan-2023.txt.gz",          
    "3LT": "9K-LT-march-2021.txt.gz",         
    "multiome": "737K-arc-v1.txt.gz",
    "flex_v1": "737K-fixed-rna-profiling.txt.gz", 
    "flex_v2": "737K-flex-v2.txt.gz"            
  }

  sc_gex_wl_dict = {
    "3v2_5v1_5v2": "cellranger_gex_barcode_whitelist_3v2_5v1_5v2.txt", 
    "3v3": "cellranger_gex_barcode_whitelist_3v3.txt", 
    "3v4": "cellranger_gex_barcode_whitelist_3v4.txt", 
    "5v3": "cellranger_gex_barcode_whitelist_5v3.txt", 
    "3LT": "cellranger_gex_barcode_whitelist_3LT.txt", 
    "multiome": "cellranger_gex_barcode_whitelist_multiome_gex.txt",
    "flex_v1": "cellranger_gex_barcode_whitelist_flex_v1.txt",
    "flex_v2": "cellranger_gex_barcode_whitelist_flex_v2.txt"
  }

  translation_cr_wl_dict = {
    "3v3": "3M-february-2018_NXT.txt.gz",          
    "3v4": "3M-3pgex-may-2023_NXT.txt.gz",         
    "3LT": "9K-LT-march-2021.txt.gz"           
  }
  
  translation_sc_wl_dict = {
    "3v3": "cellranger_whitelist_translation_3v3.txt", 
    "3v4": "cellranger_whitelist_translation_3v4.txt", 
    "3LT": "cellranger_whitelist_translation_3LT.txt"
  }

  sc_hto_wl_dict = {
    "3v3": "cellranger_hto_barcode_whitelist_3v3.txt", 
    "3v4": "cellranger_hto_barcode_whitelist_3v4.txt", 
    "3LT": "cellranger_hto_barcode_whitelist_3LT.txt",
  }
  

  print("Extracting whitelists from CellRanger")
  # extract gex whitelist and translation files from cellranger and rename them
  _extract_whitelists(tar_path, cr_wls_dict, sc_gex_wl_dict, output_dir, is_translation=False)
  _extract_whitelists(tar_path, translation_cr_wl_dict, translation_sc_wl_dict, output_dir, is_translation=True)

  # extract hto whitelists from translation files (the second column of barcodes in the translation file corresponds to hto barcodes)
  print("Extracting hto whitelists")
  _get_hto_wl_from_translation(sc_hto_wl_dict, translation_sc_wl_dict, output_dir)

  # create a lookup table for all whitelists

  chem_lu_dict = {
    "3LT": ["3LT"],
    "3v2_5v1_5v2": ["3v2", "5v1", "5v2"],
    "3v3": ["3v3"],
    "3v4": ["3v4"],
    "5v3": ["5v3"],
    "multiome": ["multiome"],
    "flex_v1": ["flex_v1"],
    "flex_v2": ["flex_v2"]
    }

  chem_ord = ["3LT", "3v2", "3v3", "3v4", "5v1", "5v2", "5v3", "multiome", "flex_v1", "flex_v2"]
    
  chem_rows = []
  reverse_map = {label: dict_key for dict_key, labels in chem_lu_dict.items() for label in labels}

  for chem in chem_ord:
    dict_key = reverse_map[chem]
    chem_rows.append({
     "chemistry": chem,
     "gex_barcodes_f": sc_gex_wl_dict.get(dict_key, ""),
     "hto_barcodes_f": sc_hto_wl_dict.get(dict_key, ""),
     "translation_f": translation_sc_wl_dict.get(dict_key, "")
    })

  chem_df = pl.from_dicts(chem_rows)
  chem_df.write_csv(whitelists_lu_f)

  # cleanup
  print("Cleaning up")
  os.remove(tar_path)
  
  return
  

def _extract_whitelists(tar_path, cr_wls_dict, sc_wl_dict, output_dir, is_translation = False): 
  # extract whitelist files from cellranger and rename them
  with tarfile.open(tar_path, "r:gz") as tar:
    for member in tar.getmembers():
      full_path = member.name
      filename  = os.path.basename(full_path)

      is_in_translation_dir = "barcodes/translation/" in full_path
      if is_translation != is_in_translation_dir:
        continue

      if "lib/python/cellranger/barcodes/" in full_path:
        for key, cr_name in cr_wls_dict.items():
          if filename == cr_name and key in sc_wl_dict:
            sc_name     = sc_wl_dict[key]
            final_path  = os.path.join(output_dir, sc_name)
                  
            f_extracted = tar.extractfile(member)
            if f_extracted:
              if filename.endswith(".gz"):
                with gzip.GzipFile(fileobj=f_extracted) as gz:
                    with open(final_path, 'wb') as f_out:
                      shutil.copyfileobj(gz, f_out)
                        
              else:
                with open(final_path, 'wb') as f_out:
                  shutil.copyfileobj(f_extracted, f_out)

  return


def _get_hto_wl_from_translation(sc_hto_wl_dict, translation_sc_wl_dict, output_dir): 
  for chem in sc_hto_wl_dict.keys(): # Changed .key to .keys()
    translation_f = os.path.join(output_dir, translation_sc_wl_dict[chem])
    hto_wl_f      = os.path.join(output_dir, sc_hto_wl_dict[chem])
    
    translation_df = pl.read_csv(translation_f, has_header=False, separator= '\t')

    translation_df.select(pl.col("column_1")).write_csv(
      hto_wl_f, include_header=False
      )

  return




def get_txome_index_parameters(config):
  # initialize
  ref_txomes  = config.get('ref_txomes', {})
  SETUP_LS    = []

  # get parameters for all specified tenx genomes
  if 'tenx' in ref_txomes:
    # get just tenx list
    tenx_ls     = ref_txomes['tenx']
    for spec_tenx in tenx_ls:
      SETUP_LS.append(_get_index_parameters_tenx(spec_tenx, TENX_NAMES))

  # get parameters for all specified custom genomes
  if 'custom' in ref_txomes:
    # get just custom list
    custom_ls     = ref_txomes['custom']
    for spec_custom in custom_ls:
      SETUP_LS.append(_get_index_parameters_custom(spec_custom))

  # check no duplicate names
  setup_names = [ s['name'] for s in SETUP_LS ]
  if not len(setup_names) == len(set(setup_names)):
    raise KeyError("Duplicated genome names are not allowed!")

  # check if reference for tutorial needs to be added
  if 'mouse_2024' not in setup_names:
    spec_tmp    = {'name': "mouse_2024"}
    setup_names.append('mouse_2024')
    SETUP_LS.append( _get_index_parameters_tenx(spec_tmp, TENX_NAMES) )

  # make dictionary instead
  SETUP_LS    = dict(zip(setup_names, SETUP_LS))

  return SETUP_LS


def _get_index_parameters_tenx(spec_tenx, TENX_NAMES):

  if spec_tenx['name'] not in TENX_NAMES:
    raise ValueError(f"name of tenx reference transcriptome {spec_tenx['name']} is incorrect")
  
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


def _get_index_parameters_custom(spec_custom):
  # check if all required entries are specified
  if not os.path.isfile(spec_custom['gtf']):
    raise FileNotFoundError(f"file {spec_custom['gtf']} specified in configfile doesn't exist")
  if not os.path.isfile(spec_custom['fasta']):
    raise FileNotFoundError(f"file {spec_custom['fasta']} specified in configfile doesn't exist")

  # fill in defaults if we need to
  if not "decoys" in spec_custom:
    spec_custom["decoys"]     = True
  if not "rrnas" in spec_custom:
    spec_custom["rrnas"]      = True

  return spec_custom


def get_probe_set_parameters(config):
  probe_sets  = config.get('probe_sets', {})
  SETUP_LS    = []

  if 'tenx' in probe_sets:
    tenx_ls   = probe_sets['tenx']
    for spec_tenx in tenx_ls:
      SETUP_LS.append(_get_probe_set_parameters_tenx(spec_tenx))

  # check no duplicate names
  setup_names = [s['name'] for s in SETUP_LS]
  if not len(setup_names) == len(set(setup_names)):
    raise KeyError("Duplicated probe set names are not allowed!")

  SETUP_LS    = dict(zip(setup_names, SETUP_LS))

  return SETUP_LS


def _get_probe_set_parameters_tenx(spec_tenx):
  if spec_tenx['name'] not in TENX_PROBE_SET_NAMES:
    raise ValueError(f"Probe set name '{spec_tenx['name']}' is not valid. "
      f"Choose from: {TENX_PROBE_SET_NAMES}")

  spec_tenx['mito_str']     = TENX_PROBE_SET_MITOS[spec_tenx['name']]
  spec_tenx['is_probe_set'] = True

  return spec_tenx


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



# function that makes simpleaf index
def set_up_txome_index(scdata_dir, txome_name, fasta_f, gtf_f, index_dir, mito_str, is_prebuilt, is_tenx, has_decoy, has_rrna, n_cores):
  if is_prebuilt:
    print(f"'is_prebuilt' is True, value is {is_prebuilt}")
  else:
    print(f"'is_prebuilt' is False, value is {is_prebuilt}")

  # create output directories
  ref_dir   = os.path.join(scdata_dir, 'reference_transcriptomes', txome_name)
  idx_dir   = os.path.join(scdata_dir, 'alevin_fry_home', 'ref_txomes', txome_name)
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


def set_up_probe_set_index(scdata_dir, probe_set_name, n_cores):
  # download probe set CSV
  probe_csv_dir = os.path.join(scdata_dir, 'probe_sets', probe_set_name)
  os.makedirs(probe_csv_dir, exist_ok=True)

  probe_csv_url = URLS_10X_PROBE_SETS[probe_set_name]
  probe_csv_f   = os.path.join(probe_csv_dir, f'{probe_set_name}_probe_set.csv')

  print(f'Downloading probe set CSV for {probe_set_name}')
  subprocess.run(['wget', '-O', probe_csv_f, probe_csv_url], check=True)

  # build alevin index from probe set
  idx_dir = os.path.join(scdata_dir, 'alevin_fry_home', 'probe_sets', probe_set_name)
  os.makedirs(idx_dir, exist_ok=True)

  _build_index_from_probe_set(probe_set_name, probe_csv_f, idx_dir, n_cores)

  # create gene info file
  gene_info_f = os.path.join(probe_csv_dir, f'{probe_set_name}_gene_info.txt.gz')
  _make_probe_set_gene_info_file(probe_set_name, probe_csv_f, gene_info_f)

  # save yaml with parameters
  mito_str  = TENX_PROBE_SET_MITOS[probe_set_name]
  yaml_f    = os.path.join(idx_dir, f'{probe_set_name}_index_params.yaml')
  _make_probe_set_index_params_yaml(yaml_f, probe_set_name, probe_csv_f, idx_dir, mito_str, gene_info_f)

  print(f'Completed making probe set index for {probe_set_name} in {scdata_dir}.')

  return


def _build_index_from_probe_set(probe_set_name, probe_csv_f, idx_dir, n_cores):
  print(f'Building alevin index from probe set {probe_set_name}')

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
  TMP_IDX_DIR="${{ALEVIN_FRY_HOME}}/{probe_set_name}"

  # simpleaf index from probe set CSV
  simpleaf index \
    --output {idx_dir} \
    --probe-csv {probe_csv_f} \
    --overwrite --threads {n_cores} \
    --work-dir ${{TMP_IDX_DIR}}

  # tidy up
  rm -rf ${{ALEVIN_FRY_HOME}}
  """

  print(bash_script)
  subprocess.run(bash_script, shell=True, executable='/bin/bash', check=True)

  return


def _make_probe_set_index_params_yaml(yaml_f, probe_set_name, probe_csv_f, idx_dir, mito_str, gene_info_f):
  param_ls = {
    "reference":    probe_set_name,
    "probe_csv_f":  probe_csv_f,
    "index_dir":    idx_dir,
    "gtf_f":        None,
    "gene_info_f":  gene_info_f,
    "mito_str":     mito_str,
    "has_decoy":    False,
    "has_rrna":     False,
    "is_prebuilt":  False,
    "is_tenx":      True,
    "is_probe_set": True
  }
  with open(yaml_f, 'w') as f:
    yaml.dump(param_ls, f)

  return


def _make_probe_set_gene_info_file(probe_set_name, probe_csv_f, gene_info_f):
  print(f'Creating gene info file for probe set {probe_set_name}')

  # read probe csv (skip comment lines), get one row per unique gene
  probe_df = pl.read_csv(probe_csv_f, comment_prefix='#')
  gene_df  = probe_df.select(['gene_id', 'gene_name']).unique()

  # genet gene annotations and chromosome info from the relevant gtf txt file
  gtf_txt_url = URLS_PROBE_GTF_TXTS[probe_set_name]
  print(f'Reading gene annotation from {gtf_txt_url}')
  gtf_txt_df  = pl.read_csv(gtf_txt_url, separator='\t')

  # join probe genes with gtf annotation on ensembl_id
  gene_info_df = (
    gene_df.join(
      gtf_txt_df.select(['ensembl_id', 'gene_type', 'chromosome', 'start', 'end', 'strand', 'width']),
      left_on='gene_id',
      right_on='ensembl_id',
      how='left'
    )
    .rename({'gene_id': 'ensembl_id', 'gene_name': 'symbol'})
    .with_columns(
      (pl.col('symbol') + '_' + pl.col('ensembl_id')).alias('gene_id')
    )
    .select(['gene_id', 'ensembl_id', 'symbol', 'gene_type', 'chromosome', 'start', 'end', 'strand', 'width'])
  )

  # save
  buf = gene_info_df.write_csv(separator='\t').encode()
  with gzip.open(gene_info_f, 'wb') as f:
    f.write(buf)

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
  TMP_IDX_DIR="${{ALEVIN_FRY_HOME}}/{ref_txome}"

  # simpleaf index
  simpleaf index \
    --output {idx_dir} \
    --fasta {fasta_f} \
    --gtf {gtf_f} \
    {decoy_flag} --overwrite --threads {n_cores} \
    --work-dir ${{TMP_IDX_DIR}} \
    --use-piscem

  # tidy up
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
    "reference":    ref_txome, 
    "fasta_f":      fasta_f, 
    "index_dir":    index_dir, 
    "gtf_f":        gtf_f, 
    "gene_info_f":  gtf_txt_f, 
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
  # check yamls exist
  assert all([ os.path.isfile(f) for f in yaml_fs]), "not all yaml files exist"

  # get data we need from each one
  df_data = []
  for yaml_f in yaml_fs:
    with open(yaml_f) as f:
      yaml_ls = yaml.load(f, Loader=yaml.FullLoader)
      df_data.append({
        "reference":   yaml_ls.get('reference') or yaml_ls.get('ref_txome'),
        "mito_str":    yaml_ls.get('mito_str'),
        "gene_info_f": yaml_ls.get('gene_info_f') or yaml_ls.get('gtf_txt_f')
      })

  # create DataFrame, save
  params_df = pl.DataFrame(df_data)
  params_df.write_csv(csv_f)


if __name__ == "__main__":
  # set up parser
  parser      = argparse.ArgumentParser()
  subparsers  = parser.add_subparsers(dest="function_name", help="Name of the function to run")

  # parser for get_scprocess_data
  getdata     = subparsers.add_parser('get_scprocess_data')
  getdata.add_argument('scdata_dir', type=str)
  getdata.add_argument('cr_url', type = str)
  getdata.add_argument('wl_lu_f', type=str)
  getdata.add_argument('cr_version_f', type=str)

  # parsers for set_up_txome_index
  get_af      = subparsers.add_parser('set_up_txome_index')
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

  # parser for set_up_probe_set_index
  get_ps      = subparsers.add_parser('set_up_probe_set_index')
  get_ps.add_argument('scdata_dir', type = str)
  get_ps.add_argument('probe_set_name', type = str, help = 'probe set name (e.g. human_v1)')
  get_ps.add_argument('cores', type = int)

  # decide which function
  args = parser.parse_args()
  if args.function_name == 'get_scprocess_data':
    get_scprocess_data(args.scdata_dir, args.cr_url, args.wl_lu_f, args.cr_version_f)
  elif args.function_name == 'set_up_txome_index':
    set_up_txome_index(args.scdata_dir, args.genome, args.fasta_f, args.gtf_f, args.index_dir, args.mito_str, 
      _safe_boolean(args.is_prebuilt), _safe_boolean(args.is_tenx), 
      _safe_boolean(args.has_decoy), _safe_boolean(args.has_rrna), args.cores)
  elif args.function_name == 'save_index_params_csv':
    save_index_params_csv(args.csv_f, args.yaml_fs)
  elif args.function_name == 'set_up_probe_set_index':
    set_up_probe_set_index(args.scdata_dir, args.probe_set_name, args.cores)
  else:
    parser.print_help()

