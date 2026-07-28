import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import yaml

from scripts.scprocess_utils import (
  check_setup_before_running_scprocess,
  get_conda_prefix,
)


class TestCondaPrefix(unittest.TestCase):
  def test_uses_scprocess_local_default(self):
    with tempfile.TemporaryDirectory() as tmp:
      self.assertEqual(
        get_conda_prefix({'user': {}}, Path(tmp)),
        Path(tmp) / '.snakemake' / 'conda',
      )

  def test_profile_setting_takes_precedence(self):
    with tempfile.TemporaryDirectory() as tmp:
      profile_dir = Path(tmp) / 'profile'
      profile_dir.mkdir()
      with (profile_dir / 'config.yaml').open('w') as stream:
        yaml.safe_dump({'conda-prefix': '/scratch/user/scprocess_envs'}, stream)

      self.assertEqual(
        get_conda_prefix({'user': {'profile_dir': profile_dir}}, Path(tmp)),
        Path('/scratch/user/scprocess_envs'),
      )

  def test_local_run_sets_cores_and_scprocess_local_conda_prefix(self):
    with tempfile.TemporaryDirectory() as tmp:
      root = Path(tmp)
      scprocess_dir = root / 'scprocess'
      scdata_dir = root / 'scdata'
      scprocess_dir.mkdir()
      scdata_dir.mkdir()
      for name in [
        'cellranger_ref', 'gmt_pathways', 'marker_genes', 'xgboost',
        'alevin_fry_home',
      ]:
        (scdata_dir / name).mkdir()
      (scdata_dir / 'index_parameters.csv').touch()
      (scdata_dir / 'scprocess_setup.yaml').write_text('user: {}\n')

      setup_cfg = {'user': {'local_cores': 3}}
      with (
        patch.dict('os.environ', {'SCPROCESS_DATA_DIR': str(scdata_dir)}),
        patch(
          'scripts.scprocess_utils.check_setup_config',
          return_value=setup_cfg,
        ),
      ):
        _, extraargs, returned_cfg = check_setup_before_running_scprocess(
          scprocess_dir, []
        )

      expected_prefix = scprocess_dir / '.snakemake' / 'conda'
      self.assertEqual(
        extraargs,
        ['--cores', '3', '--conda-prefix', str(expected_prefix)],
      )
      self.assertEqual(returned_cfg['_conda_prefix'], expected_prefix)


if __name__ == '__main__':
  unittest.main()
