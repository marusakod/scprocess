import os
import argparse
import pathlib
import subprocess

# set up
import collections
import tempfile
import gc
import random
import json
import polars as pl
import numpy as np
import yaml
import shutil

def map_fastqs_to_counts(run, af_dir, demux_type, what, af_home_dir, 
    where, R1_fs, R2_fs, threads, t2g_f, index_dir, wl_lu_f, tenx_chemistry = 'none', exp_ori = 'none', whitelist_f = 'none'):
  # make output directory, in subdirectory if multiplexed samples
  out_dir   = f"{af_dir}/af_{run}"
  if demux_type == "hto":
    if what == "rna":
      out_dir = f"{out_dir}/rna"
    elif what == "hto":
      out_dir = f"{out_dir}/hto"
  os.makedirs(out_dir, exist_ok = True)
  print('made out_dir')

  # set up simpleaf
  os.environ["ALEVIN_FRY_HOME"] = af_home_dir
  subprocess.run(["simpleaf", "set-paths"])

  # if arvados, download to temp files
  on_arvados  = not os.path.exists(where)
  if on_arvados:
    # set up tmp directory
    tmp_dir     = f"{af_dir}/.tmp_fastqs_{run}_{what}"
    prefix      = f"{run}_{what}"
    os.makedirs(tmp_dir, exist_ok = True)

    # download files from Arvados
    print('downloading files from arvados')
    arv_uuid    = where
    R1_fs       = [ _download_arvados_file_as_tempfile(arv_uuid, f, tmp_dir, prefix, i, "R1", threads) for i, f in enumerate(R1_fs) ]
    R2_fs       = [ _download_arvados_file_as_tempfile(arv_uuid, f, tmp_dir, prefix, i, "R2", threads) for i, f in enumerate(R2_fs) ]
  else:
    R1_fs       = [ os.path.join(where, f) for f in R1_fs]
    R2_fs       = [ os.path.join(where, f) for f in R2_fs]

  # get whitelist lookup file
  wl_lu_dt = pl.read_csv(wl_lu_f)

  if tenx_chemistry == 'none': 
    wl_overlap_dt = _get_whitelist_overlap(R1_fs, wl_lu_f, wl_lu_dt)
    # check for which barcode whitelist the overlap is the highest
    max_overlap = max(wl_overlap_dt['overlap'])
    if max_overlap < 0.7:
      raise Warning(f'Maximum overlap ob barcodes is {max_overlap:.1%}, 10x chemistry guess might be incorrect')
    
    sel_wl_dt = wl_overlap_dt.filter(pl.col('overlap') == max_overlap)
    whitelist_f = sel_wl_dt['barcodes_f_full'][0]
    if sel_wl_dt.height == 1:
      sample_chem = sel_wl_dt['chemistry'][0]
      tenx_chemistry = '10xv3'
      if sample_chem =='5v3':
        exp_ori = 'rc'
      else:
        exp_ori = 'fw'

      pct_mapped = "" # no mapping to downsampled data needs to be done
    else: # if selected whitelist corresponds to multiple chemistries with different orrientation, do mapping on downsampled data
      
      # make directory for temporary output files
      tmp_out_dir = f'{out_dir}/tmp_mapping'
      os.makedirs(tmp_out_dir, exist_ok=True)
      
      sub_R1_f, sub_R2_f = _subset_fastqs(tmp_out_dir, R1_fs, R2_fs)
      chem_opts = set(sel_wl_dt['chemistry'])
      if chem_opts == set(['3v2', '5v1', '5v2']):
        tenx_chemistry = '10xv2'
      else:
        tenx_chemistry = '10xv3'
      
      # map downsampled fastqs 
      _run_simpleaf_quant(tmp_out_dir, [sub_R1_f], [sub_R2_f], threads, index_dir, 
        tenx_chemistry, 'fw', t2g_f, whitelist_f)
      
      # infer read orientation
      exp_ori, pct_mapped = _infer_read_orientation(tmp_out_dir)
      
      # remove temporary mapping results and downsampled fastqs
      shutil.rmtree(tmp_out_dir)
       
      # get sample chemisty
      if exp_ori == 'fw': 
        sample_chem = '3v2'
      else:
        sample_chem = '5v1/5v2'
  
  else:
    pct_mapped  = ""
    max_overlap = ""
    # get sample chemistry based on barcode whitelist and exp_ori
    chem_opts = (wl_lu_dt
      .filter(pl.col('barcodes_f') == os.path.basename(whitelist_f))
      .get_column('chemistry').to_list())
    if set(chem_opts) == set(['3v2', '5v1', '5v2']):
      if exp_ori == 'fw': 
        sample_chem = '3v2'
      else:
        sample_chem = '5v1/5v2'
    else:
      sample_chem = chem_opts[0]
 

  # do quantification
  extra = ["--no-piscem"] if what == "hto" else []
  _run_simpleaf_quant(out_dir, R1_fs, R2_fs, threads, index_dir,
    tenx_chemistry, exp_ori, t2g_f, whitelist_f,  extra_args=extra
  )
  
  # save yaml with chemistry stats only if what is rna
  if what == 'rna':
    chem_stats_f = os.path.join(out_dir, 'chemistry_statistics.yaml')
    # get translation file
    trans_fs = (wl_lu_dt
     .filter(pl.col("barcodes_f") == os.path.basename(whitelist_f))
     .get_column("translation_f").to_list())

    if trans_fs[0] is None:
     trans_f = ""
    else:
     trans_f  = f'{os.path.dirname(whitelist_f)}/{trans_fs[0]}'

    chem_stats = {
     "run": run, 
     "selected_whitelist": whitelist_f, 
     "selected_translation_f": trans_f,
     "selected_whitelist_overlap": max_overlap, 
     "selected_ori": exp_ori, 
     "percent_mapped_fw": pct_mapped, 
     "selected_tenx_chemistry": sample_chem, 
     "selected_af_chemisty": tenx_chemistry
    }
  
    with open(chem_stats_f, "w") as f:
     yaml.safe_dump(chem_stats, f)

  # tidy up any temp fastq files
  if on_arvados:
    for f in R1_fs:
      os.unlink(f)
    for f in R2_fs:
      os.unlink(f)
    os.rmdir(tmp_dir)



def _run_simpleaf_quant(out_dir, R1_fs, R2_fs, threads, index_dir, chemistry, ori, t2g_f, wl_f, extra_args=None):
  
  simpleaf_cmd  = [
    "simpleaf", "quant", 
    "--reads1", ",".join(R1_fs), 
    "--reads2", ",".join(R2_fs),
    "--threads", f"{threads}", 
    "--index", index_dir, 
    "--chemistry", chemistry, 
    "--resolution", "cr-like", 
    "--expected-ori", ori, 
    "--t2g-map", t2g_f, 
    "--unfiltered-pl", wl_f,
    "--min-reads", "1", 
    "--output", out_dir
    ]
  if extra_args:
    simpleaf_cmd.extend(extra_args)
    
  subprocess.run(simpleaf_cmd, check=True)



def _download_arvados_file_as_tempfile(arv_uuid, f, tmp_dir, prefix, i, read, threads):

  # create a temporary file to store the data
  temp_file   = pathlib.Path(tmp_dir) / f"{prefix}.{i}.{read}.fastq.gz"

  # write the contents of the arvados file-like object to the temporary file
  print(f"  downloading {f} from arvados as tmp file {temp_file.name}")
  subprocess.run(["arv-get", f"{arv_uuid}/{f}", str(temp_file), "--threads", str(threads)])

  return str(temp_file)


def _subset_fastqs(out_dir, R1_fs, R2_fs, smpl_size = 100000):
  # check how many R1 + R2 file pairs
  if len(R1_fs) > 1:
    random.seed(12346)
    idx = random.sample(range(len(R1_fs)), 1)[0]  
    R1_f = R1_fs[idx]
    R2_f = R2_fs[idx]
  else:
    R1_f = R1_fs[0]
    R2_f = R2_fs[0]
  
  # get names for downsampled files
  R1_f_base = os.path.basename(R1_f)
  R2_f_base = os.path.basename(R2_f)
  
  # get fastq dir
  fastq_dir = os.path.dirname(R1_f)
  sub_R1_f = f'{out_dir}/downsampled_{R1_f_base}'
  sub_R2_f = f'{out_dir}/downsampled_{R2_f_base}'
  subprocess.run(["seqkit", "head", "-n", f"{smpl_size}", R1_f, "-o", sub_R1_f], check=True)
  subprocess.run(["seqkit", "head", "-n", f"{smpl_size}", R2_f, "-o", sub_R2_f], check=True)
  
  return sub_R1_f, sub_R2_f

# af_res_dir should be a directory where temporary alevin inputs for inference of chemistry are stored
def _infer_read_orientation(af_res_dir):
  
  # check the meta_info.json generated by simpleaf
  json_path = os.path.join(af_res_dir, "af_map", "map_info.json")
  with open(json_path) as f:
    data = json.load(f)
    pct_mapped = data.get("percent_mapped", 0)
    
    # high mapping rate (>50%) --> mapping with fw seems correct --> likely 3'
    if pct_mapped >= 50.0:
      ori_guess = "fw"
    # low mapping rate (<15%) --> highly likely mapping should've been with rc --> likely 5'
    elif pct_mapped <= 15.0:
      ori_guess = "rc"
    # hmm 
    else:
      ori_guess = "fw" 
      
    return ori_guess, pct_mapped


def _get_whitelist_overlap(R1_fs, wl_lu_f, wl_lu_dt, sample_size = 100000):
  # randomly pick one R1 file to extract barcodes from
  random.seed(1234)
  sel_R1_f = random.sample(R1_fs, 1)[0]
    
  # get all barcode whitelist files
  wl_dt  = wl_lu_dt.select(['chemistry', 'barcodes_f'])
  wl_fs  = wl_dt['barcodes_f'].unique().to_list()

  # get directory where whitelist files are stored
  wl_dir = os.path.abspath(os.path.dirname(wl_lu_f)) 
    
  # calculate overlap of barcodes in sel_R1_f with each whitelist 
  print(f'Extracting barcodes from {sel_R1_f}')
  spell     = f"seqkit head -n {sample_size} {sel_R1_f} | seqkit subseq -r 1:16"
  spell_res = subprocess.run(spell, shell=True, capture_output=True, text=True)
  barcodes  = set(_extract_raw_seqs_from_fq(spell_res.stdout))
  n_bcs     = len(barcodes)
  print(f'Number of unique barcodes: {n_bcs}')
  
  overlap_res = []
  for wl_f in wl_fs:
    wl_f_full = f'{wl_dir}/{wl_f}'
    with open(wl_f_full, 'r') as f: 
      wl_set = {line.strip() for line in f}
      matches = sum(1 for bc in barcodes if bc in wl_set)
      overlap_pct = matches/n_bcs if n_bcs > 0 else 0
      overlap_res.append({"barcodes_f": wl_f, "barcodes_f_full": wl_f_full, "overlap": overlap_pct})
  
  # merge overlaps with chemistries
  overlap_dt = pl.DataFrame(overlap_res)
  full_dt    = wl_dt.join(overlap_dt, on = 'barcodes_f', coalesce=True, how = 'full')

  return full_dt
    

# A FASTQ file contains sequences and their associated quality scores. Each entry is structured as
#@SEQ_ID - identifier for the sequencing read
#SEQUENCE - DNA sequence
#+ - separator line (empty or SEQ_ID)
#QUALITY - quality score string for the sequence
# this function extracts DNA sequences from each entry of a FASTQ file
def _extract_raw_seqs_from_fq(fastq_all):
  entries = fastq_all.strip().split('\n@')[0:]
  seqs = []
  for e in entries:
    ls = e.split('\n')
    if len(ls) > 1:
      seqs.append(ls[1])
  return seqs     


if __name__ == "__main__":
  # get arguments
  parser  = argparse.ArgumentParser()
  parser.add_argument("run", type=str)
  parser.add_argument("--af_dir", type=str)
  parser.add_argument("--demux_type", type=str)
  parser.add_argument("--what", default="rna", type=str, choices=["rna", "hto"])
  parser.add_argument("--af_home_dir", type=str)
  parser.add_argument("--where", type=str)
  parser.add_argument("--R1_fs", nargs="+")
  parser.add_argument("--R2_fs", nargs="+")
  parser.add_argument("--threads", default=1, type=int)
  parser.add_argument("--af_index_dir", type=str)
  parser.add_argument("--wl_lu_f", type=str)
  parser.add_argument("--tenx_chemistry", type=str, default='none')
  parser.add_argument("--exp_ori", type=str, default='none')
  parser.add_argument("--whitelist_f", type=str, default='none')

  # set up some locations
  args    = parser.parse_args()
  if args.what == 'hto':
    t2g_f     = f"{args.af_dir}/t2g_hto.tsv"
    index_dir = f"{args.af_dir}/hto_index"
  else:
    t2g_f     = f"{args.af_index_dir}/index/t2g_3col.tsv"
    index_dir = f"{args.af_index_dir}/index"

  # run
  map_fastqs_to_counts(run = args.run, af_dir = args.af_dir, demux_type = args.demux_type, what = args.what, af_home_dir = args.af_home_dir, 
    where = args.where, R1_fs=args.R1_fs, R2_fs=args.R2_fs, threads=args.threads, tenx_chemistry=args.tenx_chemistry, 
    exp_ori = args.exp_ori, wl_lu_f= args.wl_lu_f, whitelist_f= args.whitelist_f, t2g_f=t2g_f, index_dir= index_dir)


