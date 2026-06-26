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
  af_chemistry = 'none', exp_ori = 'none', whitelist_f = 'none', lib_pool_dir = ''):
  # make output directory, in subdirectory if multiplexed with HTO
  out_dir   = f"{af_dir}/{lib_pool_dir}af_{run}"
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

  # resolve FASTQ paths (download from Arvados if needed)
  R1_fs, R2_fs, on_arvados, tmp_dir = _resolve_fastq_paths(
    where, R1_fs, R2_fs, arv_instance, af_dir, lib_pool_dir, run, what, threads)

  # resolve chemistry, orientation, and whitelist
  wl_lu_dt = pl.read_csv(wl_lu_f)
  chem = _resolve_chemistry(
    R1_fs, R2_fs, af_chemistry, exp_ori, whitelist_f, wl_lu_f, wl_lu_dt,
    index_dir, t2g_f, threads, out_dir, what)

  # do quantification; HTO index has a subdirectory structure from simpleaf index
  quant_index_dir = os.path.join(index_dir, "index") if what == "hto" else index_dir
  _run_simpleaf_quant(out_dir, R1_fs, R2_fs, threads, quant_index_dir,
    chem["af_chemistry"], chem["exp_ori"], t2g_f, chem["whitelist_f"]
  )

  # save chemistry stats YAML (RNA only)
  if what == 'rna':
    _write_chemistry_stats(out_dir, run, chem, wl_lu_dt)

  # tidy up any temp fastq files
  if on_arvados:
    for f in R1_fs:
      os.unlink(f)
    for f in R2_fs:
      os.unlink(f)
    os.rmdir(tmp_dir)


def _resolve_fastq_paths(where, R1_fs, R2_fs, arv_instance, af_dir, lib_pool_dir, run, what, threads):
  """Resolve FASTQ file paths, downloading from Arvados if necessary."""
  on_arvados = not os.path.exists(where)
  tmp_dir = None
  if on_arvados:
    tmp_dir   = f"{af_dir}/{lib_pool_dir}.tmp_fastqs_{run}_{what}"
    prefix    = f"{run}_{what}"
    os.makedirs(tmp_dir, exist_ok = True)
    print('downloading files from arvados')
    arv_uuid  = where
    R1_fs     = [ _download_arvados_file_as_tempfile(arv_uuid, arv_instance, f, tmp_dir, prefix, i, "R1", threads) for i, f in enumerate(R1_fs) ]
    R2_fs     = [ _download_arvados_file_as_tempfile(arv_uuid, arv_instance, f, tmp_dir, prefix, i, "R2", threads) for i, f in enumerate(R2_fs) ]
  else:
    R1_fs     = [ os.path.join(where, f) for f in R1_fs ]
    R2_fs     = [ os.path.join(where, f) for f in R2_fs ]
  return R1_fs, R2_fs, on_arvados, tmp_dir


def _resolve_chemistry(R1_fs, R2_fs, af_chemistry, exp_ori, whitelist_f, wl_lu_f, wl_lu_dt,
  index_dir, t2g_f, threads, out_dir, what):
  """Detect or validate 10x chemistry, orientation, and whitelist.

  When af_chemistry is 'none', auto-detects via barcode overlap and optional
  orientation inference. Otherwise validates the specified chemistry against
  the R1 read length.

  Returns a dict with: af_chemistry, exp_ori, whitelist_f, sample_tenx_chemistry,
  r1_length, max_overlap, cell_counts_fw, cell_counts_rc.
  """
  if af_chemistry == 'none':
    return _autodetect_chemistry(R1_fs, R2_fs, wl_lu_f, wl_lu_dt, index_dir,
      t2g_f, threads, out_dir)
  else:
    return _validate_chemistry(R1_fs, af_chemistry, exp_ori, whitelist_f,
      wl_lu_dt, what)


def _autodetect_chemistry(R1_fs, R2_fs, wl_lu_f, wl_lu_dt, index_dir, t2g_f, threads, out_dir):
  """Auto-detect chemistry via barcode overlap, read length, and orientation inference."""

  # step 1: sample barcodes from R1 and find the best-matching whitelist
  print(' checking overlap of barcodes with different whitelists')
  sel_wl_dt, whitelist_f, max_overlap = _get_whitelist_overlap(R1_fs, wl_lu_f, wl_lu_dt)

  # step 2: check R1 read length (10xv3 requires >= 28bp; shorter falls back to 10xv2)
  r1_length = _get_r1_read_length(R1_fs)
  print(f' R1 read length: {r1_length}bp')

  # step 3: determine chemistry and orientation from the matched whitelist
  if sel_wl_dt.height == 1:
    # only one chemistry uses this whitelist — unambiguous
    sample_tenx_chemistry = sel_wl_dt['chemistry'][0]
    af_chemistry          = '10xv3'
    exp_ori               = 'rc' if sample_tenx_chemistry == '5v3' else 'fw'
    cell_counts_fw        = ""
    cell_counts_rc        = ""
  else:
    # multiple chemistries share this whitelist (e.g. 3v2/5v1/5v2 all use 737K)
    # — need to infer orientation by mapping downsampled FASTQs in both directions
    print(' guessing orientation by mapping downsampled FASTQ files')
    tmp_out_dir = f'{out_dir}/tmp_mapping'
    os.makedirs(tmp_out_dir, exist_ok=True)

    sub_R1_f, sub_R2_f    = _subset_fastqs(tmp_out_dir, R1_fs, R2_fs)
    chem_opts             = set(sel_wl_dt['chemistry'])
    af_chemistry          = '10xv2' if chem_opts == set(['3v2', '5v1', '5v2']) else '10xv3'

    # map downsampled FASTQs twice (forward and reverse complement)
    for ori in ['fw', 'rc']:
      _run_simpleaf_quant(f'{tmp_out_dir}/{ori}_mapping', [sub_R1_f], [sub_R2_f], threads, index_dir,
        af_chemistry, ori, t2g_f, whitelist_f)

    # pick the orientation that yielded more quantified cells
    exp_ori, cell_counts_fw, cell_counts_rc = _infer_read_orientation(tmp_out_dir)

    # derive the specific 10x chemistry from the inferred orientation
    if chem_opts == set(['3v2', '5v1', '5v2']):
      sample_tenx_chemistry = '3v2' if exp_ori == 'fw' else '5v1/5v2'
    else:
      sample_tenx_chemistry = '3v3' if exp_ori == 'fw' else '5v3'

    shutil.rmtree(tmp_out_dir)

  # override to 10xv2 if R1 is too short for 10xv3
  if af_chemistry == '10xv3' and r1_length < 28:
    warnings.warn(f'R1 read length is {r1_length}bp (< 28bp); using 10xv2 chemistry with v3 whitelist')
    af_chemistry          = '10xv2'
    sample_tenx_chemistry = '3v2' if exp_ori == 'fw' else '5v1/5v2'

  return {
    "af_chemistry": af_chemistry, "exp_ori": exp_ori, "whitelist_f": whitelist_f,
    "sample_tenx_chemistry": sample_tenx_chemistry, "r1_length": r1_length,
    "max_overlap": max_overlap, "cell_counts_fw": cell_counts_fw, "cell_counts_rc": cell_counts_rc
  }


def _validate_chemistry(R1_fs, af_chemistry, exp_ori, whitelist_f, wl_lu_dt, what):
  """Validate a user-specified chemistry against the R1 read length."""
  f_prefix = 'gex' if what == 'rna' else 'hto'
  chem_opts = (wl_lu_dt
    .filter(pl.col(f'{f_prefix}_barcodes_f') == os.path.basename(whitelist_f))
    .get_column('chemistry').to_list())
  if set(chem_opts) == set(['3v2', '5v1', '5v2']):
    sample_tenx_chemistry = '3v2' if exp_ori == "fw" else '5v1/5v2'
  else:
    sample_tenx_chemistry = chem_opts[0]

  r1_length = _get_r1_read_length(R1_fs)
  print(f' R1 read length: {r1_length}bp')
  if af_chemistry == '10xv3' and r1_length < 28:
    warnings.warn(f'R1 read length is {r1_length}bp (< 28bp); using 10xv2 chemistry with v3 whitelist')
    af_chemistry = '10xv2'
    sample_tenx_chemistry = '3v2' if exp_ori == 'fw' else '5v1/5v2'

  return {
    "af_chemistry": af_chemistry, "exp_ori": exp_ori, "whitelist_f": whitelist_f,
    "sample_tenx_chemistry": sample_tenx_chemistry, "r1_length": r1_length,
    "max_overlap": "", "cell_counts_fw": "", "cell_counts_rc": ""
  }


def _write_chemistry_stats(out_dir, run, chem, wl_lu_dt):
  """Write chemistry_statistics.yaml with detection/validation results."""
  chem_stats_f = os.path.join(out_dir, 'chemistry_statistics.yaml')

  # look up translation file and HTO whitelist from the GEX whitelist
  trans_fs = (wl_lu_dt
    .filter(pl.col("gex_barcodes_f") == os.path.basename(chem["whitelist_f"]))
    .get_column("translation_f").to_list())
  hto_whitelist_fs = (wl_lu_dt
    .filter(pl.col("gex_barcodes_f") == os.path.basename(chem["whitelist_f"]))
    .get_column("hto_barcodes_f").to_list())

  if trans_fs[0] is None:
    trans_f         = ""
    hto_whitelist_f = ""
  else:
    wl_dir          = os.path.dirname(chem["whitelist_f"])
    trans_f         = f'{wl_dir}/{trans_fs[0]}'
    hto_whitelist_f = f'{wl_dir}/{hto_whitelist_fs[0]}'

  chem_stats = {
    "run": run,
    "r1_read_length": chem["r1_length"],
    "selected_gex_whitelist": chem["whitelist_f"],
    "selected_hto_whitelist": hto_whitelist_f,
    "selected_translation_f": trans_f,
    "selected_whitelist_overlap": chem["max_overlap"],
    "selected_ori": chem["exp_ori"],
    "n_cells_fw": chem["cell_counts_fw"],
    "n_cells_rc": chem["cell_counts_rc"],
    "selected_tenx_chemistry": chem["sample_tenx_chemistry"],
    "selected_af_chemistry": chem["af_chemistry"]
  }

  with open(chem_stats_f, "w") as f:
    yaml.safe_dump(chem_stats, f)


def map_flex_fastqs_to_counts(run, af_dir, af_home_dir, where,
  R1_fs, R2_fs, threads, index_dir, af_chemistry, whitelist_f, probset_f, probe_bc_f,
  arv_instance=None, lib_pool_dir='', geometry=None):
  out_dir   = f"{af_dir}/{lib_pool_dir}af_{run}/flex"
  os.makedirs(out_dir, exist_ok = True)
  print('made out_dir')

  # set up simpleaf
  os.environ["ALEVIN_FRY_HOME"] = af_home_dir
  subprocess.run(["simpleaf", "set-paths"])

  # if arvados, download to temp files
  on_arvados  = not os.path.exists(where)
  if on_arvados:
    # set up tmp directory
    tmp_dir     = f"{af_dir}/.tmp_fastqs_{run}"
    prefix      = f"{run}"
    os.makedirs(tmp_dir, exist_ok = True)

    # download files from Arvados
    print('downloading files from arvados')
    arv_uuid    = where
    R1_fs       = [ _download_arvados_file_as_tempfile(arv_uuid, arv_instance, f, tmp_dir, prefix, i, "R1", threads) for i, f in enumerate(R1_fs) ]
    R2_fs       = [ _download_arvados_file_as_tempfile(arv_uuid, arv_instance, f, tmp_dir, prefix, i, "R2", threads) for i, f in enumerate(R2_fs) ]
  else:
    R1_fs       = [ os.path.join(where, f) for f in R1_fs]
    R2_fs       = [ os.path.join(where, f) for f in R2_fs]

  # run simpleaf multiplex-quant
  simpleaf_cmd = [
    "simpleaf", "multiplex-quant", 
    "--reads1", ",".join(R1_fs), 
    "--reads2", ",".join(R2_fs),
    "--threads", f"{threads}", 
    "--index", index_dir,
    "--probe-set", probset_f,
    "--chemistry", af_chemistry,
    "--cell-bc-list", whitelist_f,
    "--sample-bc-list", probe_bc_f,
    "--usa", 
    "--expected-ori", "fw",
    "--min-reads", "1",
    "--output", out_dir
  ]

  if geometry:
    simpleaf_cmd.extend(["--geometry", geometry])

  subprocess.run(simpleaf_cmd, check=True)

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
  """Sample barcodes from R1 reads and find the best-matching whitelist.

  Returns (sel_wl_dt, whitelist_f, max_overlap): the filtered rows for the
  best whitelist (may have multiple rows if several chemistries share it),
  the full path to the whitelist file, and the overlap fraction.
  """
  # randomly pick one R1 file to extract barcodes from
  random.seed(1234)
  sel_R1_f  = random.sample(R1_fs, 1)[0]

  # get all barcode whitelist files
  wl_dt     = wl_lu_dt.select(['chemistry', 'gex_barcodes_f'])
  wl_fs     = wl_dt['gex_barcodes_f'].unique().to_list()

  # get directory where whitelist files are stored
  wl_dir    = os.path.abspath(os.path.dirname(wl_lu_f))

  # sample barcodes (first 16bp of R1) and compute overlap with each whitelist
  print(f'Extracting barcodes from {sel_R1_f}')
  spell     = f"seqkit head -n {sample_size} {sel_R1_f} | seqkit subseq -r 1:16"
  spell_res = subprocess.run(spell, shell=True, capture_output=True, text=True)
  barcodes  = set(_extract_raw_seqs_from_fq(spell_res.stdout))
  n_bcs     = len(barcodes)
  if n_bcs == 0:
    raise RuntimeError(f"No barcodes extracted from {sel_R1_f}. The file may be empty or corrupted.")
  print(f'Number of unique barcodes: {n_bcs}')

  # compute fraction of sampled barcodes found in each whitelist
  overlap_res = []
  for wl_f in wl_fs:
    wl_f_full = f'{wl_dir}/{wl_f}'
    with open(wl_f_full, 'r') as f:
      wl_set      = {line.strip() for line in f}
      matches     = sum(1 for bc in barcodes if bc in wl_set)
      overlap_pct = matches/n_bcs
      overlap_res.append({"gex_barcodes_f": wl_f, "gex_barcodes_f_full": wl_f_full, "overlap": overlap_pct})

  # merge overlaps with chemistries and pick the best match
  overlap_dt  = pl.DataFrame(overlap_res)
  full_dt     = wl_dt.join(overlap_dt, on = 'gex_barcodes_f', coalesce=True, how = 'full')
  max_overlap = max(full_dt['overlap'])
  if max_overlap < 0.7:
    warnings.warn(f'Maximum overlap of barcodes is {max_overlap:.1%}, 10x chemistry guess might be incorrect')
  sel_wl_dt   = full_dt.filter(pl.col('overlap') == max_overlap)
  whitelist_f = sel_wl_dt['gex_barcodes_f_full'][0]

  return sel_wl_dt, whitelist_f, max_overlap
    

def _extract_raw_seqs_from_fq(fastq_all):
  """Extract DNA sequences from FASTQ-formatted text.

  Each FASTQ entry has four lines:
    @SEQ_ID  - identifier for the sequencing read
    SEQUENCE - DNA sequence
    +        - separator line
    QUALITY  - quality score string

  This function splits on record boundaries and returns the SEQUENCE line
  from each entry.
  """
  entries = fastq_all.strip().split('\n@')[0:]
  seqs = []
  for e in entries:
    ls = e.split('\n')
    if len(ls) > 1:
      seqs.append(ls[1])
  return seqs


def _get_r1_read_length(R1_fs):
  random.seed(1234)
  sel_R1_f  = random.sample(R1_fs, 1)[0]
  spell     = f"seqkit head -n 1 {sel_R1_f} | seqkit seq -s"
  spell_res = subprocess.run(spell, shell=True, capture_output=True, text=True)
  seq       = spell_res.stdout.strip()
  if not seq:
    raise RuntimeError(f"Could not read R1 sequence from {sel_R1_f}")
  return len(seq)


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  subparsers = parser.add_subparsers(dest='command', required=True)

  # map polya fastqs to counts
  p_map = subparsers.add_parser('map_fastqs_to_counts')
  p_map.add_argument("run", type=str)
  p_map.add_argument("--af_dir", type=str)
  p_map.add_argument("--demux_type", type=str, default="none")
  p_map.add_argument("--what", default="rna", type=str, choices=["rna", "hto"])
  p_map.add_argument("--af_home_dir", type=str)
  p_map.add_argument("--where", type=str)
  p_map.add_argument("--R1_fs", nargs="+")
  p_map.add_argument("--R2_fs", nargs="+")
  p_map.add_argument("--threads", default=1, type=int)
  p_map.add_argument("--af_index_dir", type=str)
  p_map.add_argument("--wl_lu_f", type=str)
  p_map.add_argument("--af_chemistry", type=str, default='none')
  p_map.add_argument("--exp_ori", type=str, default='none')
  p_map.add_argument("--whitelist_f", type=str, default='none')
  p_map.add_argument("--arv_instance", type=str, default=None)
  p_map.add_argument("--lib_pool_dir", type=str, default='')

  # map flex fastqs to counts
  p_flex = subparsers.add_parser('map_flex_fastqs_to_counts')
  p_flex.add_argument("run", type=str)
  p_flex.add_argument("--af_dir", type=str)
  p_flex.add_argument("--af_home_dir", type=str)
  p_flex.add_argument("--where", type=str)
  p_flex.add_argument("--R1_fs", nargs="+")
  p_flex.add_argument("--R2_fs", nargs="+")
  p_flex.add_argument("--threads", default=1, type=int)
  p_flex.add_argument("--af_index_dir", type=str)
  p_flex.add_argument("--af_chemistry", type=str)
  p_flex.add_argument("--gex_whitelist_f", type=str)
  p_flex.add_argument("--probeset_f", type=str)
  p_flex.add_argument("--probe_bc_f", type=str)
  p_flex.add_argument("--arv_instance", type=str, default=None)
  p_flex.add_argument("--geometry", type=str, default=None)
  p_flex.add_argument("--lib_pool_dir", type=str, default='')

  args = parser.parse_args()

  if args.command == 'map_fastqs_to_counts':
    if args.what == 'hto':
      t2g_f     = f"{args.af_dir}/t2g_hto.tsv"
      index_dir = f"{args.af_dir}/hto_index"
    else:
      t2g_f     = f"{args.af_index_dir}/index/t2g_3col.tsv"
      index_dir = f"{args.af_index_dir}/index"
      
    map_fastqs_to_counts(run=args.run, af_dir=args.af_dir, demux_type=args.demux_type,
      what=args.what, af_home_dir=args.af_home_dir, where=args.where,
      R1_fs=args.R1_fs, R2_fs=args.R2_fs, threads=args.threads, af_chemistry=args.af_chemistry,
      exp_ori=args.exp_ori, wl_lu_f=args.wl_lu_f, whitelist_f=args.whitelist_f, t2g_f=t2g_f,
      index_dir=index_dir, arv_instance=args.arv_instance, lib_pool_dir=args.lib_pool_dir)

  elif args.command == 'map_flex_fastqs_to_counts':
    index_dir = f"{args.af_index_dir}/index"
    map_flex_fastqs_to_counts(run=args.run, af_dir=args.af_dir, af_home_dir=args.af_home_dir,
      where=args.where, R1_fs=args.R1_fs, R2_fs=args.R2_fs, threads=args.threads,
      index_dir=index_dir, af_chemistry=args.af_chemistry, whitelist_f=args.gex_whitelist_f,
      probset_f=args.probeset_f, probe_bc_f=args.probe_bc_f, arv_instance=args.arv_instance,
      lib_pool_dir=args.lib_pool_dir, geometry=args.geometry)



