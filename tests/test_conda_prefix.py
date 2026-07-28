import tempfile
import unittest
from pathlib import Path

import yaml

from scripts.scprocess_utils import get_conda_prefix


class TestCondaPrefix(unittest.TestCase):
  def test_uses_shared_scprocess_data_default(self):
    with tempfile.TemporaryDirectory() as tmp:
      self.assertEqual(
        get_conda_prefix({'user': {}}, Path(tmp)),
        Path(tmp) / 'conda_envs',
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


if __name__ == '__main__':
  unittest.main()
