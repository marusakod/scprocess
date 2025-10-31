import os
import argparse
import pathlib
import subprocess

# set up
import collections
import tempfile
import gc

def map_fastqs_to_counts(run, af_dir, demux_type, what, af_home_dir, 
    where, R1_fs, R2_fs, threads, af_index_dir, tenx_chemistry, 
    exp_ori, whitelist_f, t2g_f, index_dir):
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

  # do quantification
  simpleaf_cmd  = [
    "simpleaf", "quant", 
    "--reads1", ",".join(R1_fs), 
    "--reads2", ",".join(R2_fs),
    "--threads", f"{threads}", 
    "--index", index_dir, 
    "--chemistry", tenx_chemistry, 
    "--resolution", "cr-like", 
    "--expected-ori", exp_ori, 
    "--t2g-map", t2g_f, 
    "--unfiltered-pl", whitelist_f,
    "--min-reads", "1", 
    "--output", out_dir
    ]
  if what == "hto":
    simpleaf_cmd.append("--no-piscem")
  subprocess.run(simpleaf_cmd)

  # tidy up any temp fastq files
  if on_arvados:
    for f in R1_fs:
      os.unlink(f)
    for f in R2_fs:
      os.unlink(f)
    os.rmdir(tmp_dir)


def _download_arvados_file_as_tempfile(arv_uuid, f, tmp_dir, prefix, i, read, threads):
  # create a temporary file to store the data
  temp_file   = pathlib.Path(tmp_dir) / f"{prefix}.{i}.{read}.fastq.gz"

  # write the contents of the arvados file-like object to the temporary file
  print(f"  downloading {f} from arvados as tmp file {temp_file.name}")
  subprocess.run(["arv-get", f"{arv_uuid}/{f}", str(temp_file), "--threads", str(threads)])

  return str(temp_file)


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
  parser.add_argument("--tenx_chemistry", type=str)
  parser.add_argument("--exp_ori", type=str)
  parser.add_argument("--whitelist_f", type=str)

  # set up some locations
  args    = parser.parse_args()
  if args.what == 'hto':
    t2g_f     = f"{args.af_dir}/t2g_hto.tsv"
    index_dir = f"{args.af_dir}/hto_index"
  else:
    t2g_f     = f"{args.af_index_dir}/index/t2g_3col.tsv"
    index_dir = f"{args.af_index_dir}/index"

  # run
  map_fastqs_to_counts(args.run, args.af_dir, args.demux_type, args.what, args.af_home_dir, 
    args.where, args.R1_fs, args.R2_fs, args.threads, args.af_index_dir, args.tenx_chemistry, 
    args.exp_ori, args.whitelist_f, t2g_f, index_dir)
