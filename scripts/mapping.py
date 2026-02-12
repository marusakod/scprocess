import argparse
import pathlib
import subprocess
import os

# set up
import collections
import tempfile
import random
import json
import polars as pl
import numpy as np
import yaml
import shutil
import warnings

# import relevant libraries for arvados download
import arvados
from concurrent.futures import ThreadPoolExecutor, as_completed
from queue import Queue
import threading


def map_fastqs_to_counts(run, af_dir, demux_type, what, af_home_dir, where, 
  R1_fs, R2_fs, threads, t2g_f, index_dir, wl_lu_f, arv_instance = None,
  tenx_chemistry = 'none', exp_ori = 'none', whitelist_f = 'none'):
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
    R1_fs       = [ _download_arvados_file_as_tempfile(arv_uuid, arv_instance, f, tmp_dir, prefix, i, "R1", threads) for i, f in enumerate(R1_fs) ]
    R2_fs       = [ _download_arvados_file_as_tempfile(arv_uuid, arv_instance, f, tmp_dir, prefix, i, "R2", threads) for i, f in enumerate(R2_fs) ]
  else:
    R1_fs       = [ os.path.join(where, f) for f in R1_fs]
    R2_fs       = [ os.path.join(where, f) for f in R2_fs]

  # get whitelist lookup file
  wl_lu_dt = pl.read_csv(wl_lu_f)

  if tenx_chemistry == 'none': 
    print(' checking overlap of barcodes with different whitelists')
    wl_overlap_dt = _get_whitelist_overlap(R1_fs, wl_lu_f, wl_lu_dt)
    # check for which barcode whitelist the overlap is the highest
    max_overlap = max(wl_overlap_dt['overlap'])
    if max_overlap < 0.7:
      warnings.warn(f'Maximum overlap ob barcodes is {max_overlap:.1%}, 10x chemistry guess might be incorrect')
    
    sel_wl_dt = wl_overlap_dt.filter(pl.col('overlap') == max_overlap)
    whitelist_f = sel_wl_dt['barcodes_f_full'][0]
    
    if sel_wl_dt.height == 1:
      sample_chem = sel_wl_dt['chemistry'][0]
      tenx_chemistry = '10xv3'
      exp_ori = 'rc' if sample_chem =='5v3' else 'fw'
      
      cell_counts_fw = ""
      cell_counts_rc = "" # no mapping to downsampled data needs to be done
    else: # if selected whitelist corresponds to multiple chemistries with different orrientation, do mapping on downsampled data
      print(' guessing orientation by mapping downsampled FASTQ files')

      # make directory for temporary output files
      tmp_out_dir = f'{out_dir}/tmp_mapping'
      os.makedirs(tmp_out_dir, exist_ok=True)
      
      sub_R1_f, sub_R2_f = _subset_fastqs(tmp_out_dir, R1_fs, R2_fs)
      chem_opts = set(sel_wl_dt['chemistry'])
      tenx_chemistry = '10xv2' if chem_opts == set(['3v2', '5v1', '5v2']) else '10xv3'
      
      # map downsampled fastqs 2x
      for ori in ['fw', 'rc']:
        _run_simpleaf_quant(f'{tmp_out_dir}/{ori}_mapping', [sub_R1_f], [sub_R2_f], threads, index_dir, 
          tenx_chemistry, ori , t2g_f, whitelist_f)
        
      # infer read orientation
      exp_ori, cell_counts_fw, cell_counts_rc = _infer_read_orientation(tmp_out_dir)

      # remove temporary mapping results and downsampled fastqs
      shutil.rmtree(tmp_out_dir)
       
      # get sample chemisty
      sample_chem = '3v2' if exp_ori == 'fw' else '5v1/5v2'
  
  else:
    cell_counts_fw = ""
    cell_counts_rc = ""
    max_overlap    = ""
    # get sample chemistry based on barcode whitelist and exp_ori
    chem_opts = (wl_lu_dt
      .filter(pl.col('barcodes_f') == os.path.basename(whitelist_f))
      .get_column('chemistry').to_list())
    if set(chem_opts) == set(['3v2', '5v1', '5v2']):
      sample_chem = '3v2' if exp_ori == "fw" else '5v1/5v2'
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
      "n_cells_fw": cell_counts_fw, 
      "n_cells_rc": cell_counts_rc, 
      "selected_tenx_chemistry": sample_chem, 
      "selected_af_chemistry": tenx_chemistry
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


def _download_arvados_file_as_tempfile(arv_uuid, arv_instance, f, tmp_dir, prefix, i, read, threads):
  """
  Download a file from Arvados using the Python API with multithreading support.
  
  Args:
    arv_uuid: Arvados collection UUID
    f: File name within the collection
    tmp_dir: Temporary directory for storing the file
    prefix: Prefix for the temporary file
    i: Index for the file
    read: Read type (e.g., "R1" or "R2")
    threads: Number of threads to use for concurrent buffering and I/O
  
  Returns:
    Path to the downloaded temporary file as a string
  """

  # Create temporary file path
  temp_file = pathlib.Path(tmp_dir) / f"{prefix}.{i}.{read}.fastq.gz"

  print(f"  downloading {f} from arvados as tmp file {temp_file.name}")

  try:
    # set up arvados access
    arv_token   = os.environ["ARVADOS_API_TOKEN"]
    arv_client  = arvados.api('v1', host = f'api.{arv_instance}.roche.com',
      token = arv_token, insecure = True, num_retries = 2 )
    
    # Open the collection
    collection = arvados.collection.Collection(arv_uuid, api_client=arv_client)
    
    # Download file using multithreaded copy
    _download_file_multithreaded(collection, f, str(temp_file), threads)
    
    return str(temp_file)
  
  except Exception as e:
    raise RuntimeError(f"Failed to download {f} from Arvados collection {arv_uuid}: {str(e)}")


def _download_file_multithreaded(collection, src_path, dest_path, num_threads, chunk_size=1024*1024):
  """
  Download a file from an Arvados collection to disk using multithreading.
  
  Uses a queue-based approach where reader threads pull chunks from the source
  and writer threads push them to disk. This allows for efficient buffering and
  parallelization of I/O operations.
  
  Args:
    collection: Arvados collection object
    src_path: Path to file within the collection
    dest_path: Destination file path on disk
    num_threads: Number of threads to use (min 2: 1 reader, 1 writer)
    chunk_size: Size of chunks to read/write (default 1 MB)
  """
  # Ensure at least 2 threads (1 for reading, 1 for writing)
  num_threads = max(2, num_threads)
  
  chunk_queue = Queue(maxsize=num_threads * 2)
  exceptions = []
  
  def reader():
    """Read chunks from source file in the collection."""
    try:
      with collection.open(src_path, 'rb') as src:
        while True:
          chunk = src.read(chunk_size)
          if not chunk:
            chunk_queue.put(None)  # Signal end of file
            break
          chunk_queue.put(chunk)
    except Exception as e:
      chunk_queue.put(('ERROR', str(e)))
      exceptions.append(e)
  
  def writer():
    """Write chunks to destination file."""
    try:
      with open(dest_path, 'wb') as dst:
        while True:
          item = chunk_queue.get()
          if item is None:
            break
          if isinstance(item, tuple) and item[0] == 'ERROR':
            raise RuntimeError(f"Read error: {item[1]}")
          dst.write(item)
    except Exception as e:
      exceptions.append(e)
      raise
  
  # Start reader and writer threads
  reader_thread = threading.Thread(target=reader, daemon=False)
  writer_thread = threading.Thread(target=writer, daemon=False)
  
  reader_thread.start()
  writer_thread.start()
  
  # Wait for both threads to complete
  reader_thread.join()
  writer_thread.join()
  
  # Check for any exceptions
  if exceptions:
    raise RuntimeError(f"Error during multithreaded download: {exceptions[0]}")


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


def _infer_read_orientation(af_res_dir):
  
  # get quant.json files from simpleaf
  json_paths = {
    'fw': os.path.join(af_res_dir, 'fw_mapping', 'af_quant', 'quant.json'),
    'rc': os.path.join(af_res_dir, 'rc_mapping', 'af_quant', 'quant.json')
  }
    
  cell_counts = {}
  
  for ori, path in json_paths.items():
    if not os.path.exists(path):
      raise FileNotFoundError(f"Mapping file not found: {path}")
            
    with open(path) as f:
      data = json.load(f)
      cell_counts[ori] = data.get("num_quantified_cells", 0)

  if cell_counts['fw'] > cell_counts['rc']:
    ori_guess = "fw"
  elif cell_counts['rc'] > cell_counts['fw']:
    ori_guess = "rc"
  else:
    raise ValueError(f"Ambiguous orientation: same number of cells quantified with both orientations")
    
  return ori_guess, cell_counts['fw'], cell_counts["rc"]


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
  parser.add_argument("--arv_instance", type=str, default=None)

  # set up some locations
  args    = parser.parse_args()
  if args.what == 'hto':
    t2g_f     = f"{args.af_dir}/t2g_hto.tsv"
    index_dir = f"{args.af_dir}/hto_index"
  else:
    t2g_f     = f"{args.af_index_dir}/index/t2g_3col.tsv"
    index_dir = f"{args.af_index_dir}/index"

  # run
  map_fastqs_to_counts(run = args.run, af_dir = args.af_dir, demux_type = args.demux_type, 
    what = args.what, af_home_dir = args.af_home_dir, where = args.where, 
    R1_fs=args.R1_fs, R2_fs=args.R2_fs, threads=args.threads, tenx_chemistry=args.tenx_chemistry, 
    exp_ori = args.exp_ori, wl_lu_f= args.wl_lu_f, whitelist_f= args.whitelist_f, t2g_f=t2g_f, 
    index_dir= index_dir, arv_instance = args.arv_instance)


