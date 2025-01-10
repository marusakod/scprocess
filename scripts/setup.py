# functions for setting up scprocess_data

import pandas as pd
import pyranges as pr
import gzip
import argparse
import os
import re
import glob
import subprocess
import numpy as np

#configfile = '/Users/marusa/Projects/scprocess_test/configs/config-setup-template.yaml'

def get_genome_params(GENOME_NAMES, FASTA_FS, GTF_FS, INDEX_DIRS, MITO_STRS, DECOYS, RRNAS, SCPROCESS_DATA_DIR):

  # make dictionaries with all parameters that will be updated later
  fasta_dict    = dict(zip(GENOME_NAMES, FASTA_FS))
  gtf_dict      = dict(zip(GENOME_NAMES, GTF_FS))
  idx_dict      = dict(zip(GENOME_NAMES, INDEX_DIRS))
  mito_str_dict = dict(zip(GENOME_NAMES, MITO_STRS))
  decoys_dict   = dict(zip(GENOME_NAMES, DECOYS))
  rrnas_dict    = dict(zip(GENOME_NAMES, RRNAS))
  
  # crete directory for reference genomes inside scprocess data
  ref_dir = os.path.join(SCPROCESS_DATA_DIR, 'reference_genomes')
  os.makedirs(ref_dir, exist_ok = True)

  #create directory for alevin indices inside scprocess data
  af_dir = os.path.join(SCPROCESS_DATA_DIR, 'alevin_fry_home')
  os.makedirs(af_dir, exist_ok= True)

  # sort out predefined genomes
  pre_gnomes = {'human_2020', 'human_2024','mouse_2020', 'mouse_2024'}
  all_gnomes = set(GENOME_NAMES)

# which predefined genomes to download
  names = list(pre_gnomes & all_gnomes)
  
  if(len(names) != 0):

    refs_10x_links_dict ={
      'human_2020': "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz", 
      'human_2024': "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz", 
      'mouse_2020': "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz", 
      'mouse_2024': "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
    }

    # download all required genomes from 10x
    for n in names:
        # if rna and decoy download a prebuild index
        dcoys = decoys_dict[n]
        rrnas = rrnas_dict[n]
      
        if dcoys and rrnas:
          fasta_dict[n] = None
          gtf_dict[n] = None

        elif not dcoys and rrnas:
          gnome_fs = build_10x_genome_w_rrna(ref_dir = ref_dir, name = n)
          fasta_dict[n] = gnome_fs[1]
          gtf_dict[n] = gnome_fs[0]

        else:
          l = refs_10x_links_dict[n]
          gnome_fs = get_predefined_gnomes(ref_dir= ref_dir, name = n, link = l )
          # add fasta_f and gtf_f to dicts
          fasta_dict[n] =  gnome_fs[1]
          gtf_dict[n] = gnome_fs[0]
    
    # define mito strings for predefined genomes
    for  n in names:
      mito_str = "^MT-" if 'human_' in n else "^mt-"
      mito_str_dict[n] = mito_str

# sort out custom genomes
  cust_names = list(all_gnomes - pre_gnomes)
  if(len(cust_names) != 0):
    print('Sorting out locations for custom genomes')
    # create directories for custom genomes in ref_dir and symlinks for fasta and gtf files 
    for cn in cust_names:
      # get original fasta
      cn_fasta = fasta_dict[cn]
      # get original gtf
      cn_gtf = gtf_dict[cn]
      # get index dir
      cn_idx_dir = idx_dict[cn]

      cn_ref_dir = os.path.join(ref_dir, cn)
      os.makedirs(cn_ref_dir, exist_ok = True)
      gtf_sym = os.path.join(cn_ref_dir, 'genes.gtf')
     
      if not os.path.exists(gtf_sym):
        os.symlink(cn_gtf, gtf_sym)
      
      gtf_dict[cn] = gtf_sym
      
      # create symlink to index if exists othewise symlink to fasta
      if cn_idx_dir is not None:
        print(f"Using premade index for {cn}. Ignoring 'decoys' parameter")
        cn_af_idx_dir = os.path.join(af_dir, n)
        os.makedirs(cn_af_idx_dir, exist_ok=True)

        create_idx_symlinks(cn_idx_dir, cn_af_idx_dir)
        
      else:
        fasta_sym = os.path.join(cn_ref_dir, 'genes.fa')
        fasta_dict[cn] = fasta_sym
        if not os.path.exists(fasta_sym):
          os.symlink(cn_fasta, fasta_sym)

        fasta_dict[cn] = fasta_sym
  

  # make a txt file from all gtf files (or download from zenodo)
  
  gtf_txt_dict = {}

  for n, gtf in gtf_dict.items():

    if gtf is None:
       
      zenodo_gtf_urls ={ 
        'human_2020': "https://zenodo.org/records/14247195/files/human_2020_rrna_genes_gtf.txt.gz", 
        'human_2024': "https://zenodo.org/records/14247195/files/human_2024_rrna_genes_gtf.txt.gz", 
        'mouse_2020': "https://zenodo.org/records/14247195/files/mouse_2020_rrna_genes_gtf.txt.gz", 
        'mouse_2024': "https://zenodo.org/records/14247195/files/mouse_2024_rrna_genes_gtf.txt.gz"
      }

      url = zenodo_gtf_urls[n]

      print(f'Downloading gtf txt file for {n}')
      out_dir = os.path.join(ref_dir, n)
      subprocess.run(f"wget -P {out_dir} {url}", shell=True)
      # rename gtf txt file
      txt_f     = glob.glob( out_dir + '/*genes_gtf.txt.gz')

      assert len(txt_f) == 1, \
       f"More than one gtf.txt file found in {out_dir}"
      
      out_txt_f = os.path.join(out_dir, 'genes_gtf.txt.gz')
      os.rename(txt_f[0], out_txt_f)
      gtf_txt_dict.update({n:out_txt_f})

    else:
      print(f'Creating txt files from gtf for {n}')
      # define output file
      out_txt_f = os.path.join(os.path.dirname(gtf), 'genes_gtf.txt.gz')
      # convert .gtf to .txt
      save_gtf_as_txt(gtf, out_txt_f)
      # update dictionary 
      gtf_txt_dict.update({n:out_txt_f})


  # create one dataframe from all dictionaries 
  PARAMS_DF = pd.DataFrame([fasta_dict, gtf_dict, gtf_txt_dict, mito_str_dict, decoys_dict, rrnas_dict])
  PARAMS_DF = PARAMS_DF.T
  PARAMS_DF.columns = ['fasta_f', 'gtf_f', 'gtf_txt_f', 'mito_str', 'decoy', 'rrna']
  PARAMS_DF.index.name = 'genome_name'
  PARAMS_DF.reset_index(inplace=True)
  
  out_csv_f = os.path.join(SCPROCESS_DATA_DIR, 'setup_parameters.csv')
  PARAMS_DF.to_csv(out_csv_f, index = False)
  
  print(f'Completed writing genomes in {SCPROCESS_DATA_DIR}')
  return 
  


def build_10x_genome_w_rrna(ref_dir, name):
   
   print(f"Creating {name} 10x genome with rRNAs")
 
   # get bash script to download gtf and fasta
   bash_f = f"scripts/build_10x_style_genomes/build_10x_style_{name}_genome_w_rRNAs.sh"
   bash_f = os.path.realpath(bash_f)

   # create a new directory for the specified genome and switch to it
   gnome_dir = os.path.join(ref_dir, name)
   os.makedirs(gnome_dir, exist_ok=True)
   os.chdir(gnome_dir)

   # run bash script that downloads gtf and 
   subprocess.run(bash_f, shell=True)

   print('Done!')

   return [os.path.join(gnome_dir, 'genes.gtf'), os.path.join(gnome_dir, 'genome.fa')]



# get refrence genomes from 10x
def get_predefined_gnomes(ref_dir, name, link):

  print(f'Downloading {name} genome from 10x')
  
  # get tarball 
  tball = link.rsplit('/', 1)[-1]
  
  # create a new directory for the specified genome and switch to it
  gnome_dir = os.path.join(ref_dir, name)
  os.makedirs(gnome_dir, exist_ok = True)
  os.chdir(gnome_dir)
        
  # download reference into that dir
  subprocess.run(f'wget {link}', shell=True)
        
  # list tar contents
  tar_out= subprocess.run(f'tar -tzf {tball}', shell=True, capture_output=True, text=True)
  tball_all_fs = tar_out.stdout.splitlines()
        
  # look for genome.fa and genes.gtf (or genes.gtf.gz) files and extract the
  patt = re.compile(r'(genome\.fa|genes\.gtf(\.gz)?)$') 
  tball_filt_fs = [file for file in tball_all_fs if patt.search(file)]

  for t in tball_filt_fs:
    ext_spell = ['tar', '--extract', '--file=' + tball, '--strip-components=2', t]
    subprocess.run(ext_spell, check=True)

  # unzip
  gnome_fs = os.listdir()
  for f in gnome_fs:
    if f.endswith('.gtf.gz'):
      subprocess.run(f'gunzip {f}', shell=True, capture_output=False)

  # delete the tarball
  os.remove(tball)
  print('Done!')

  return [os.path.join(gnome_dir, 'genes.gtf'), os.path.join(gnome_dir, 'genome.fa')]




# download  marker genes, xgboost objects, barcode whitelists and gmt pathways from github
def get_scprocess_data(scprocess_data_dir):

  print('Downloading data from scprocessData github repo')
  # switch to scprocess data dir
  os.chdir(scprocess_data_dir)

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



# convert genes.gtf file to txt file
def save_gtf_as_txt(gtf_f, gtf_txt_f):
  
    gtf = pr.read_gtf(gtf_f, as_df = True)

    # restrict to just genes
    gene_annots = gtf[gtf['Feature'] == 'gene']
    # remove some columnns
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

    # calculate width and add gene_id column
    gene_annots['width'] = gene_annots['end'] - gene_annots['start'] + 1
    gene_annots['gene_id'] = gene_annots['symbol'] + '_' + gene_annots['ensembl_id']
    
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


def create_idx_symlinks(src_dir, target_dir):
    for root, _, files in os.walk(src_dir): 
        # get relative path from src_dir to current directory
        rel_path = os.path.relpath(root, src_dir)
        
        # map relative path to target_dir
        target_root = os.path.join(target_dir, rel_path)
        
        os.makedirs(target_root, exist_ok=True)
        
        # create symlinks for files
        for f in files:
            src_file_path = os.path.join(root, f)
            target_file_path = os.path.join(target_root, f)
            
            if not os.path.exists(target_file_path):
                os.symlink(src_file_path, target_file_path)
            else:
                print(f"Symlink already exists: {target_file_path}")



def list_of_strings(arg):
    return arg.split(',')

def list_of_bools(arg):
    str_ls = arg.split(',')
    bool_list = [True if s == 'True' else False for s in str_ls]
    return bool_list


  # make some functions executable from the command line (make_af_idx, get_genome_params, get_scprocess_data)


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  subparsers = parser.add_subparsers(dest="function_name", help="Name of the function to run")

    # parser for get_scprocess_data
  parser_getdata = subparsers.add_parser('get_scprocess_data')
  parser_getdata.add_argument('scp_data_dir', type=str)

    # parsers for get_genome_params
  parser_getgnomes = subparsers.add_parser('get_genome_params')
  parser_getgnomes.add_argument('genomes', type=list_of_strings, help='list with all genome names')
  parser_getgnomes.add_argument('fasta_fs', type=list_of_strings, help='list with paths to all fasta files')
  parser_getgnomes.add_argument('gtf_fs', type=list_of_strings, help='list with paths to all gtf files')
  parser_getgnomes.add_argument('idx_dirs', type =list_of_strings, help='list with paths to prebuild alevin indices' )
  parser_getgnomes.add_argument('mito_str', type=list_of_strings, help='list with all mitochondrial gene identifiers')
  parser_getgnomes.add_argument('decoys', type=list_of_bools, help='list with bool values definiing whether or not to use decoys when building indices with simpleaf')
  parser_getgnomes.add_argument('rrnas', type=list_of_bools, help='list with bool values definiing whether or not to include ribosomal rrnas for simpleaf index')
  parser_getgnomes.add_argument('scp_data_dir', type=str)

  args = parser.parse_args()

  if args.function_name == 'get_scprocess_data':
      get_scprocess_data(args.scp_data_dir)
  elif args.function_name == 'get_genome_params':
      get_genome_params(args.genomes, args.fasta_fs, args.gtf_fs, args.idx_dirs, args.mito_str, args.decoys, args.rrnas, args.scp_data_dir)
  else:
    parser.print_help()

