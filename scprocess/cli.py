#!/usr/bin/env python3

# import utils
import argparse
import datetime
import filecmp
import os
import shlex
import shutil
import subprocess
import textwrap
from pathlib import Path

import polars as pl
import yaml

from scprocess import utils as scprocess_utils


def _resolve_conda_prefix(setup_cfg: dict) -> Path:
  user = setup_cfg.get("user", {})
  if "conda_prefix" in user:
    return Path(user["conda_prefix"])
  base = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
  return base / "scprocess"


def _clean_conda_env():
  """Return a copy of os.environ without few conda related variables
  to minimize conda related clashes.
  """
  env = os.environ.copy()
  env.pop("CONDA_EXE", None)

  # Remove exported conda bash function from the environment
  for key in list(env.keys()):
    if key.startswith("BASH_FUNC_conda"):
      del env[key]

  return env


# scprocess setup
def run_setup(scprocess_dir, snakefile, ranger_url, dryrun, extraargs):
  print("doing some checks on inputs")

  # get scprocess_data_dir; should be defined in .bashrc
  scdata_dir = Path(os.getenv("SCPROCESS_DATA_DIR"))
  if scdata_dir:
    print(f"SCPROCESS_DATA_DIR is set to: {scdata_dir}")
  else:
    raise ValueError("SCPROCESS_DATA_DIR is not defined in .bashrc")

  # check that scdata_dir is a directory
  if not os.path.isdir(scdata_dir):
    raise FileNotFoundError("SCPROCESS_DATA_DIR is not a directory")

  # check that config file exists
  setup_f = scdata_dir / "scprocess_setup.yaml"
  if not os.path.exists(setup_f):
    raise FileNotFoundError(f"Config file {setup_f} does not exist")

  # get validated config
  with open(setup_f, "r") as f:
    setup_cfg = yaml.safe_load(f)
  schema_f = Path(sc_dir) / "resources/schemas/setup.schema.json"
  if not schema_f.is_file():
    raise FileNotFoundError("setup schema file not found")
  setup_cfg = scprocess_utils.check_setup_config(setup_cfg, schema_f, scprocess_dir)

  log_dir = scdata_dir / ".log/scprocess"
  (log_f, log_header) = _get_main_log(setup_f, log_dir, dryrun)

  with open(log_f, "w") as f:
    f.write(log_header)

  # handle slurm logs
  if "profile_dir" in setup_cfg["user"]:
    extraargs.append("--workflow-profile")
    extraargs.append(str(setup_cfg["user"]["profile_dir"]))

    profile_f = setup_cfg["user"]["profile_dir"] / "config.yaml"
    with open(profile_f, "r") as f:
      profile = yaml.safe_load(f)

    if profile.get("executor") == "slurm":
      cluster_log_dir = scdata_dir / ".log/cluster"
      os.makedirs(cluster_log_dir, exist_ok=True)
      extraargs.append("--slurm-logdir")
      extraargs.append(str(cluster_log_dir))
  else:
    extraargs.append("--cores")
    extraargs.append(str(setup_cfg["user"]["local_cores"]))

  # get cellranger version
  if ranger_url == "":
    # setup was done already, look for version file
    ranger_version_f = scdata_dir / "cellranger_ref/.cellranger_version"
    if not os.path.exists(ranger_version_f):
      raise ValueError("Please provide a valid CellRanger download link")

    with open(ranger_version_f, "r") as f:
      ranger_version = f.read().strip()
  else:
    ranger_version = scprocess_utils.check_ranger_url(ranger_url)

    # set up commands
  cmd_list = [
    "snakemake",
    "--directory",
    str(scdata_dir),
    "--snakefile",
    str(snakefile),
    "--configfile",
    str(setup_f.resolve()),
    "--config",
    f"cellranger_url={ranger_url}",
    f"cellranger_version={ranger_version}",
    "--rerun-triggers",
    "mtime",
    "params",
    "--rerun-incomplete",
    "--show-failed-logs",
    "--printshellcmds",
    "--use-conda",
    "--conda-prefix",
    str(_resolve_conda_prefix(setup_cfg)),
  ]
  cmd_list.extend(extraargs)

  # run
  cmd_str = " ".join(map(shlex.quote, map(str, cmd_list)))
  strip_colors = "sed 's/\\x1b\\[[0-9;]*[a-zA-Z]//g'"
  final_cmd = f"script -q -e -c {shlex.quote(cmd_str)} /dev/null | tee /dev/stderr | {strip_colors} >> {log_f}"

  print("starting scprocess setup")
  print(f"snakemake command:\n{cmd_str}")
  print(f"Logging progress to: {log_f}")
  subprocess.run(final_cmd, check=True, shell=True, env=_clean_conda_env())


# scprocess newproj
def new_proj(sc_dir, name, where_path, create_subdirs, create_config):
  # define template path
  src_dir = sc_dir / "resources" / "newproj_files"
  where_path = Path(where_path)
  proj_dir = where_path / name

  # check if where_path is a valid path
  if not where_path.exists():
    raise ValueError(f"Error: The specified directory '{where_path}' does not exist.")
  if not where_path.is_dir():
    raise ValueError(f"Error: '{where_path}' is not a directory.")

  # check if project already exists
  if proj_dir.exists():
    raise ValueError("project already exists")

  try:
    # create project directory
    os.mkdir(proj_dir)
    os.chdir(proj_dir)

    # create required directories
    dir_ls = [".log", "data", "code", "analysis", "output", "public"]
    for d in dir_ls:
      os.mkdir(d)

    # create additional subdirectories if requested
    if create_subdirs:
      os.mkdir(os.path.join(proj_dir, "data", "fastqs"))
      os.mkdir(os.path.join(proj_dir, "data", "metadata"))

    # define required files and locations
    f_ls = [
      "_workflowr.yml",
      ".gitignore",
      ".gitattributes",
      "_site.yml",
      "custom.css",
      "about.Rmd",
      "index.Rmd",
      "license.Rmd",
      ".nojekyll",
    ]
    loc_ls = [
      ".",
      ".",
      ".",
      "analysis",
      "analysis",
      "analysis",
      "analysis",
      "analysis",
      "public",
    ]
    assert len(f_ls) == len(loc_ls)

    # copy template files
    for f, l in zip(f_ls, loc_ls):
      subprocess.call(["cp", os.path.join(src_dir, f), os.path.join(proj_dir, l, f)])

    # Modify _site.yml
    site_f = os.path.join(proj_dir, "analysis", "_site.yml")
    with open(site_f, "r") as file:
      site_yml_txt = file.read()
    site_yml_txt = site_yml_txt.replace("proj_template", name)
    with open(site_f, "w") as file:
      file.write(site_yml_txt)

    # copy project file
    proj_ext = ".Rproj"
    subprocess.call(
      [
        "cp",
        os.path.join(src_dir, "proj_template" + proj_ext),
        os.path.join(proj_dir, name + proj_ext),
      ]
    )

    # create .Rprofile file
    rprofile_txt = """## This makes sure that R loads the workflowr package
    ## automatically, every time the project is loaded
    if (requireNamespace("workflowr", quietly = TRUE)) {
    message("Loading .Rprofile for the current workflowr project")
    library("workflowr")
    } else {
    message("workflowr package not installed, please run install.packages(\'workflowr\') to use the workflowr functions")
    }
    """
    rprofile_txt = textwrap.dedent(rprofile_txt)
    rprofile_f = os.path.join(proj_dir, ".Rprofile")
    with open(rprofile_f, "w") as file:
      file.write(rprofile_txt)

    # create config file if requested
    if create_config is not None:
      # import something we need
      from datetime import date

      # set up files etc
      config_f = os.path.join(proj_dir, f"config-{name}.yaml")
      today = date.today()
      date_stamp = today.strftime("%Y-%m-%d")
      fastq_dir = "data/fastqs" if create_subdirs else ""
      meta_dir = "data/metadata" if create_subdirs else ""

      # check that config file exists
      scdata_dir = os.getenv("SCPROCESS_DATA_DIR")
      setup_f = os.path.join(scdata_dir, "scprocess_setup.yaml")
      if not os.path.exists(setup_f):
        raise FileNotFoundError(f"Config file {setup_f} does not exist")

      # get validated config
      with open(setup_f, "r") as f:
        setup_cfg = yaml.safe_load(f)
      schema_f = Path(sc_dir) / "resources/schemas/setup.schema.json"
      if not schema_f.is_file():
        raise FileNotFoundError("setup schema file not found")
      setup_cfg = scprocess_utils.check_setup_config(setup_cfg, schema_f, sc_dir)

      # get name and affiliation from setup config if they are there
      your_name = setup_cfg["user"].get("your_name", "")
      affiliation = setup_cfg["user"].get("affiliation", "")
      arv_instance = setup_cfg["arvados"].get("arv_instance", "")
      int_use_gpu = setup_cfg["user"].get("int_use_gpu", True)

      # define config text
      proj_txt = f"""
        project:
          proj_dir: {proj_dir}
          fastq_dir: {fastq_dir}
          full_tag: {name}
          short_tag:
          your_name: {your_name}
          affiliation: {affiliation}
          sample_metadata: {meta_dir}/
          ref_txome:
          date_stamp: "{date_stamp}"
        """
      if arv_instance != "":
        arv_txt = f"""  arv_instance: {arv_instance}
        """
      else:
        arv_txt = ""

      if not int_use_gpu:
        integration_txt = f"""integration: 
          int_use_gpu: {str(int_use_gpu).lower()}
        """
      else:
        integration_txt = ""

      if "sc" in create_config:
        qc_txt = """qc:
          qc_max_mito: 0.1
          qc_min_splice: 0.10
          qc_max_splice: 0.99
        """
      if "sn" in create_config:
        qc_txt = """qc:
          qc_max_mito: 0.1
          qc_max_splice: 0.75
        """
      if "multiplex" in create_config:
        multi_txt = """multiplexing:
          demux_type:
        """
      else:
        multi_txt = ""

      # join together
      config_txt = proj_txt + arv_txt + integration_txt + multi_txt + qc_txt
      config_txt = textwrap.dedent(config_txt).strip()

      # write to file
      with open(config_f, "w") as file:
        file.write(config_txt)

  except Exception as e:
    # if an error happens, remove the directory we started making
    if os.path.exists(proj_dir):
      print(f"Error encountered. Cleaning up directory: {proj_dir}")
      shutil.rmtree(proj_dir)
      raise e  # Re-raise the error so the user knows what went wrong
  return


# scprocess run
def run_scprocess(configfile, snakefile, rule, extraargs, doindex, dryrun):
  print("doing some checks on inputs")

  # housekeeping
  config_path = Path(configfile).resolve()
  scprocess_dir = Path(__file__).parent

  # do some checks
  (scdata_dir, extraargs, setup_cfg) = (
    scprocess_utils.check_setup_before_running_scprocess(scprocess_dir, extraargs)
  )
  lm_f = scprocess_dir / "resources/snakemake/resources_lm_params_2025-12-16.csv"

  # get validated config
  with open(config_path, "r") as f:
    config = yaml.safe_load(f)
  schema_f = scprocess_dir / "resources/schemas/config.schema.json"
  if not schema_f.is_file():
    raise FileNotFoundError("schema file not found")
  config = scprocess_utils.check_config(config, schema_f, scdata_dir, scprocess_dir)

  # define master log file
  proj_dir = Path(config["project"]["proj_dir"])
  log_dir = proj_dir / ".log/scprocess"
  (log_f, log_header) = _get_main_log(config_path, log_dir, dryrun)
  # Write the header to the file immediately
  with open(log_f, "w") as f:
    f.write(log_header)

  # change cluster logs directory if executor is slurm
  if "profile_dir" in setup_cfg["user"]:
    profile_f = setup_cfg["user"]["profile_dir"] / "config.yaml"
    with open(profile_f, "r") as f:
      profile = yaml.safe_load(f)

    if profile["executor"] == "slurm":
      cluster_log_dir = proj_dir / ".log/cluster"

      extraargs.append("--slurm-logdir")
      extraargs.append(str(cluster_log_dir))

  # get lists of parameters
  (RUN_PARAMS, _) = scprocess_utils.get_run_parameters(config, scdata_dir)
  RUNS = list(RUN_PARAMS.keys())
  (BATCH_PARAMS, BATCH_VAR, SAMPLES) = scprocess_utils.get_batch_parameters(
    config, RUNS, scdata_dir
  )
  BATCHES = list(BATCH_PARAMS.keys())
  _ = scprocess_utils.prep_resource_params(config, schema_f, lm_f, RUN_PARAMS, BATCHES)
  _, _ = scprocess_utils.get_runs_to_batches(config, RUNS, BATCHES, BATCH_VAR)
  _ = scprocess_utils.get_labeller_parameters(config, schema_f, scdata_dir)

  # check config is ok for the rule we want to run
  scprocess_utils.check_config_ok_for_rule(config, rule)

  # print which files will be processed
  print("starting scprocess")
  print(f"\nprocessing {len(BATCHES)} samples:")
  print("  " + ", ".join(BATCHES))

  # assemble what we'll do
  cmd_list = [
    "snakemake",
    "--directory",
    str(proj_dir),
    "--snakefile",
    str(snakefile),
    "--configfile",
    str(config_path),
    "--config",
    f"scprocess_dir={str(scprocess_dir)}",
    "--rerun-triggers",
    "mtime",
    "params",
    "--rerun-incomplete",
    "--show-failed-logs",
    "--printshellcmds",
    "--software-deployment-method",
    "conda",
    "--conda-prefix",
    str(_resolve_conda_prefix(setup_cfg)),
    "--use-apptainer",
    "--apptainer-args",
    f"--cleanenv --nv --bind /tmp,{config['project']['proj_dir']}",
    *extraargs,
    rule,
  ]

  cmd_str = " ".join(map(shlex.quote, map(str, cmd_list)))
  strip_colors = "sed 's/\\x1b\\[[0-9;]*[a-zA-Z]//g'"
  final_cmd = (
    f"script -q -e -c {shlex.quote(cmd_str)} /dev/null"
    f" | tee /dev/stderr | {strip_colors} >> {log_f}"
  )

  print(f"snakemake command:\n{cmd_str}")
  print(f"Logging progress to: {log_f}")
  subprocess.run(final_cmd, check=True, shell=True, env=_clean_conda_env())

  # render index if requested
  if doindex:
    _render_index(proj_dir, config, config_path)

  return


def _render_index(proj_dir: Path, config, config_path):
  sc_dir = Path(__file__).parent
  # 1. Locate the template and script
  template_f = (sc_dir / "resources/rmd_templates/index.Rmd.template").resolve()
  render_script = sc_dir / "scripts/render_htmls.R"

  # unpack from config - proj_dir already passed as parameter
  rmd_dir = proj_dir / "analysis"
  docs_dir = proj_dir / "public"

  # 2. Collect config file paths for inclusion in index
  config_files = {"config_f": str(config_path)}

  # Add custom_sample_params if it exists
  if "custom_sample_params" in config["project"]:
    config_files["custom_sample_params_f"] = str(
      config["project"]["custom_sample_params"]
    )

  # Add zoom spec files if they exist
  if "zoom" in config and config["zoom"]:
    proj_dir_for_paths = Path(config["project"]["proj_dir"])
    zoom_specs_list = []
    for f in config["zoom"]:
      p = Path(f)
      if not p.is_absolute():
        p = proj_dir_for_paths / p
      zoom_specs_list.append(str(p.resolve()))
    config_files["zoom_specs"] = ",".join(zoom_specs_list)

  # 3. Build the R snippet (escaping quotes where necessary)
  # Escape backslashes in file paths for R
  config_files_escaped = {k: v.replace("\\", "\\\\") for k, v in config_files.items()}

  # Get show_arv_uuids setting (default True)
  show_arv_uuids = config["project"].get("show_arv_uuids", True)
  show_arv_uuids_r = "TRUE" if show_arv_uuids else "FALSE"

  r_code = f"""
  source('{render_script}');
  render_html(
    rule_name   = 'index',
    proj_dir    = '{proj_dir}',
    your_name   = '{config["project"]["your_name"]}',
    affiliation = '{config["project"]["affiliation"]}',
    docs_dir    = '{docs_dir}',
    short_tag   = '{config["project"]["short_tag"]}',
    full_tag    = '{config["project"]["full_tag"]}',
    date_stamp  = '{config["project"]["date_stamp"]}',
    mkr_sel_res = '{config["marker_genes"]["mkr_sel_res"]}',
    temp_f      = '{template_f}',
    rmd_f       = '{rmd_dir}/index.Rmd',
    config_f    = '{config_files_escaped["config_f"]}',
    show_arv_uuids = {show_arv_uuids_r}"""

  # Add optional config files to R code if they exist
  if "custom_sample_params_f" in config_files_escaped:
    r_code += f",\n    custom_sample_params_f = '{config_files_escaped['custom_sample_params_f']}'"

  if "zoom_specs" in config_files_escaped:
    r_code += f",\n    zoom_specs = '{config_files_escaped['zoom_specs']}'"

  r_code += """
  )
  """

  # get env
  rlibs_f = sc_dir / "envs/rlibs.yaml"
  rlibs_pins = sorted(rlibs_f.parent.glob(f"{rlibs_f.stem}.*.pin.txt"))
  env_path = _find_env_path_from_yaml(proj_dir, rlibs_f, rlibs_pins)
  if not env_path:
    raise RuntimeError(
      "Could not find a Snakemake conda environment matching rlibs.yaml "
      "or any rlibs.<platform>.pin.txt. "
      "Has the pipeline run successfully at least once?"
    )

  # actually do it
  print(f"rendering html index")
  cmd = ["conda", "run", "--prefix", env_path, "Rscript", "--vanilla", "-e", r_code]
  return subprocess.run(
    cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
  )


def _get_main_log(config_path, log_dir, dryrun):

  # get project directory and time stamp
  os.makedirs(log_dir, exist_ok=True)
  timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

  # get log file name
  log_f = "scprocess_dry_run.log" if dryrun else f"scprocess_run_{timestamp}.log"
  log_full_f = log_dir / log_f

  # get git tag
  log_header = (
    f"\n{'=' * 60}\n"
    f"RUN START:    {timestamp}\n"
    f"GIT VERSION:  {_get_version()}\n"
    f"CONFIG:       {config_path}\n"
    f"{'=' * 60}\n\n\n\n"
  )

  return (log_full_f, log_header)


def _find_env_path_from_yaml(proj_dir: Path, source_yaml, source_pin=None):
  conda_dir = proj_dir / ".snakemake" / "conda"
  source_yaml = source_yaml.resolve()

  if not os.path.exists(conda_dir):
    return None

  # Match Snakemake's copied spec files against the source YAML and optional pin file.
  # If a match is found, the environment prefix is the same basename without extension.
  source_specs = [(".yaml", source_yaml)]
  source_pins = []
  if source_pin:
    if isinstance(source_pin, (list, tuple, set)):
      source_pins = [Path(p).resolve() for p in source_pin if Path(p).exists()]
    else:
      pin_path = Path(source_pin)
      if pin_path.exists():
        source_pins = [pin_path.resolve()]
  for pin_file in source_pins:
    source_specs.append((".pin.txt", pin_file))

  for f in os.listdir(conda_dir):
    potential_match = conda_dir / f
    for suffix, source_file in source_specs:
      if f.endswith(suffix) and filecmp.cmp(
        source_file, potential_match, shallow=False
      ):
        env_path = conda_dir / f[: -len(suffix)]
        if env_path.is_dir():
          return env_path

  return None


# scprocess plotknee
def plot_interactive_knee(config_f, knee_f, sample):
  # check we can find the file
  if args.kneefile:
    knee_f = Path(args.kneefile)
    if not knee_f.is_file():
      raise FileNotFoundError(
        f"Specified knee file doesn't exist. Check for typos, but",
        " if the file doesn't exist you can generate the file by running the rule 'mapping'",
      )
  else:
    # open config file and get project directory
    with open(config_f) as f:
      config = yaml.load(f, Loader=yaml.FullLoader)

    # get required variables
    proj_dir = Path(config["project"]["proj_dir"])
    short_tag = config["project"]["short_tag"]
    date_stamp = config["project"]["date_stamp"]

    # use them to define a knee file to use
    knee_f = (
      proj_dir
      / f"output/{short_tag}_mapping/af_{sample}"
      / f"knee_plot_data_{sample}_{date_stamp}.csv.gz"
    )
    if not knee_f.is_file():
      knee_f = (
        proj_dir
        / f"output/{short_tag}_mapping/af_{sample}/rna"
        / f"knee_plot_data_{sample}_{date_stamp}.csv.gz"
      )
    if not knee_f.is_file():
      raise FileNotFoundError(
        f"Knee file for {sample} doesn't exist. You can generate it ",
        "by running the rule 'mapping'",
      )

  print("importing modules for plotting knees")
  import plotly.express as px

  # read knee file
  knee_df = pl.read_csv(knee_f, schema_overrides={"rank": pl.Float32})

  # make plot
  fig = px.scatter(
    knee_df.to_pandas(),
    x="rank",
    y="total",
    title=f"sample_id: {sample}",
    log_x=True,
    log_y=True,
    labels={"rank": "barcode rank", "total": "library size"},
  )
  fig.update_traces(marker=dict(color="black", size=5))

  # update the layout to set background colors
  fig.update_layout(
    plot_bgcolor="white",
    paper_bgcolor="white",
    xaxis=dict(gridcolor="#e5e7e9", zerolinecolor="#e5e7e9"),
    yaxis=dict(gridcolor="#e5e7e9", zerolinecolor="#e5e7e9"),
    autosize=False,
    width=1200,
    height=700,
  )

  # save plot as html file
  knee_dir = os.path.dirname(knee_f)
  plot_f = os.path.join(knee_dir, f"scprocess_knee_plot_{sample}.html")
  fig.write_html(plot_f)
  # plot(fig, auto_open=True)

  # tell user where knee plot was saved
  print(f"interactive knee plot saved here:\n  {plot_f}")


def _get_version():
  try:
    from importlib.metadata import version

    return version("scprocess")
  except Exception:
    try:
      return (
        subprocess.check_output(
          ["git", "describe", "--tags", "--always", "--dirty"],
          stderr=subprocess.DEVNULL,
        )
        .decode()
        .strip()
      )
    except Exception:
      return "unknown"


def show_version_and_hpc_info(sc_dir: Path):
  version = _get_version()
  print(f"\nscprocess version: {version}")

  # Get HPC information from scprocess_setup.yaml if available
  scdata_dir = os.getenv("SCPROCESS_DATA_DIR")
  if not scdata_dir:
    print("\nNo HPC configuration found (SCPROCESS_DATA_DIR not set)")
    return
  setup_f = Path(scdata_dir) / "scprocess_setup.yaml"
  if not setup_f.exists():
    print(f"\nNo HPC configuration found ({setup_f} does not exist)")
    return

  # see if we can read the config and display some of the details in a user-friendly way
  try:
    with open(setup_f, "r") as f:
      setup_cfg = yaml.safe_load(f)

    # Display user information
    if setup_cfg.get("user"):
      user_info = setup_cfg["user"]
      print("")
      if user_info.get("your_name"):
        print(f"  User: {user_info['your_name']}")
      if user_info.get("affiliation"):
        print(f"  Affiliation: {user_info['affiliation']}")

      if user_info.get("profile"):
        print("\nHPC Configuration:")
        profile_name = user_info["profile"]
        print(f"  Profile: {profile_name}")
        # Try to display profile config details
        profile_f = sc_dir / "profiles" / profile_name / "config.yaml"
        if profile_f.exists():
          with open(profile_f, "r") as f:
            profile_cfg = yaml.safe_load(f)
          if profile_cfg:
            print(f"  Location: {profile_f}")
            if profile_cfg.get("executor"):
              print(f"  Executor: {profile_cfg['executor']}")

    else:
      print("  No user information configured")

    print("\nReference genomes:")

    # Display reference genomes
    if setup_cfg.get("ref_txomes"):
      ref_txomes = setup_cfg["ref_txomes"]
      if ref_txomes.get("tenx"):
        print(f"  10x:")
        for ref_txomes in ref_txomes["tenx"]:
          print(f"    - {ref_txomes['name']}")
      if ref_txomes.get("custom"):
        print(f"  custom:")
        for ref_txomes in ref_txomes["custom"]:
          print(f"    - {ref_txomes['name']}")

  except Exception as e:
    print(f"\nWarning: Could not read HPC configuration: {e}")


def main():
  # define arguments
  parser = argparse.ArgumentParser(
    description="snakemake workflows for processing single cell RNAseq data."
  )
  subparsers = parser.add_subparsers(
    dest="subcommand", help="scprocess subcommands:", required=False
  )

  # subparser for the 'setup' subcommand
  setup_prsr = subparsers.add_parser("setup", help="do setup for scprocess")

  # add arguments for setup
  setup_prsr.add_argument(
    "-c",
    "--rangerurl",
    type=str,
    default="",
    help="""
      "Valid download link for CellRanger v9.0.0 or higher. Download links can be found at https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions).
      """,
  )
  setup_prsr.add_argument(
    "-n",
    "--dry-run",
    action="store_true",
    help="""
      "Dry run" execution, i.e. snakemake will print out what it would do for the setup step, but not actually do it.
      """,
  )
  setup_prsr.add_argument(
    "-E",
    "--extraargs",
    action="store",
    nargs=1,
    type=str,
    help="""
      Extra snakamake arguments. Must be provided within quotes and use an equals sign.
      For example, to have a dryrun: -E="-n"
      """,
  )

  # subparser for the 'plotknee' subcommand
  newprj_prsr = subparsers.add_parser(
    "newproj", help="create new project folder with structure that scprocess expects"
  )

  # add arguments for new project
  newprj_prsr.add_argument("name", type=str, help="Name of the project")
  newprj_prsr.add_argument(
    "-w",
    "--where",
    type=str,
    default=os.getcwd(),
    help="Where to create the project (default: current directory)",
  )
  newprj_prsr.add_argument(
    "-s",
    "--sub",
    action="store_true",
    help="Create data/fastqs and data/metadata subdirectories",
  )
  newprj_prsr.add_argument(
    "-c",
    "--config",
    type=str,
    choices=["sc", "sn", "multiplex"],
    nargs="+",
    help="Create a blank config.yml file with required scprocess parameters. Users must specify at least one option from: sc (single cell); sn (single nuclei); multiplex.",
  )

  # subparser for the 'run' subcommand
  run_prsr = subparsers.add_parser("run", help="run scprocess")

  # add arguments for scprocess
  run_prsr.add_argument(
    "configfile",
    type=str,
    help="Required. YAML file specifying what you want to run, and any non-default parameters",
  )
  run_prsr.add_argument(
    "-r",
    "--rule",
    type=str,
    default="all",
    choices=[
      "all",
      "mapping",
      "ambient",
      "demux",
      "qc",
      "hvg",
      "integration",
      "marker_genes",
      "label_celltypes",
      "zoom",
    ],
    help="""
      Leave empty to run the whole workflow, or alternatively select the rule you want to run.
      """,
  )
  run_prsr.add_argument(
    "-n",
    "--dry-run",
    action="store_true",
    help="""
      "Dry run" execution, i.e. snakemake will print out what it would do, but not actually do it.
      """,
  )
  run_prsr.add_argument(
    "-E",
    "--extraargs",
    action="store",
    nargs=1,
    type=str,
    help="""
      Extra snakamake arguments. Must be provided within quotes and given after an equals sign. 
      For example, to reduce outputs from snakemake: -E="--quiet"
      """,
  )
  run_prsr.add_argument(
    "--unlock",
    action="store_true",
    help="""
      When an scprocess run is stopped before it is finished, the snakemake directory
      may be locked. If you get an error message saying that it is locked, use this
      option to unlock it, then run scprocess as normal.
      """,
  )
  run_prsr.add_argument(
    "--create-envs",
    action="store_true",
    help="""
      Only create the conda environments needed for the workflow, without
      running any rules.
      """,
  )
  run_prsr.add_argument("--noindex", action="store_true", help="")

  # subparser for the 'plotknee' subcommand
  knee_parser = subparsers.add_parser(
    "plotknee",
    help="save interactive knee plot for a specified sample, to help specify custom parameters",
  )

  # define arguments
  knee_parser.add_argument("sample", type=str, help="Sample to be plotted.")
  group = knee_parser.add_mutually_exclusive_group()
  group.add_argument(
    "-c",
    "--configfile",
    type=str,
    nargs="?",
    default=None,
    help="Path to configuration file used for running scprocess. ",
  )
  group.add_argument(
    "-k", "--kneefile", type=str, default=None, help="Path to knee file"
  )

  # define scprocess directory
  sc_dir = Path(os.path.dirname(os.path.realpath(__file__)))

  # Parse the arguments
  args = parser.parse_args()

  # If no subcommand provided, show version and HPC info
  if args.subcommand is None:
    show_version_and_hpc_info(sc_dir)

  # create new project folder
  elif args.subcommand == "setup":
    # select snakefile corresponding to workflow
    snakefile = sc_dir / "rules/setup.smk"

    # sort out extra arguments
    extraargs = []
    if args.extraargs:
      extraargs.extend(shlex.split(args.extraargs[0]))
    if args.dry_run:
      extraargs.append("-np")

    # call bsub thing
    run_setup(sc_dir, snakefile, args.rangerurl, args.dry_run, extraargs)

  elif args.subcommand == "newproj":
    # make new project directory
    if not args.config is None:
      # check for presence of sc or sn
      has_sc = "sc" in args.config
      has_sn = "sn" in args.config

      if not (has_sc or has_sn):
        raise ValueError(
          "Error: When using --config/-c, you must specify either 'sc' (single cell) or 'sn' (single nucleus)."
        )

      if has_sc and has_sn:
        raise ValueError("Error: Please choose either 'sc' or 'sn', not both.")

    where_abs = os.path.abspath(args.where)
    new_proj(sc_dir, args.name, where_abs, args.sub, args.config)

  # run scprocess
  elif args.subcommand == "run":
    # set default rule to all
    if args.rule == "":
      rule = "all"
    else:
      rule = args.rule

    # select snakefile corresponding to workflow
    if rule == "zoom":
      snakefile = sc_dir / "rules/zoom.smk"
    else:
      snakefile = sc_dir / "rules/scprocess.smk"

    # sort out extra arguments
    extraargs = []
    noindex = args.noindex
    doindex = (
      (not noindex)
      and (not args.dry_run)
      and (not args.unlock)
      and (not args.create_envs)
    )
    if args.unlock:
      extraargs.append("--unlock")
    else:
      if args.create_envs:
        extraargs.append("--conda-create-envs-only")
      if args.extraargs:
        extraargs.extend(shlex.split(args.extraargs[0]))
      if args.dry_run:
        extraargs.append("-np")

    # call bsub thing
    run_scprocess(args.configfile, snakefile, rule, extraargs, doindex, args.dry_run)

  # plot knee for selected sample
  elif args.subcommand == "plotknee":
    # call function
    plot_interactive_knee(args.configfile, args.kneefile, args.sample)

  else:
    # if no arguments are provided at all, print the help message
    parser.print_help()


if __name__ == "__main__":
  main()
