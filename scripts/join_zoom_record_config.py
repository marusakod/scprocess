"""Record requested joined-zoom settings, independently of analytical rules.

This snapshot is not a completion flag; consult rule logs for run status.
"""
import json
from pathlib import Path

import yaml


def record_config(context, config, output):
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, 'w') as fh:
        yaml.safe_dump({'config': config, 'zoom_context': context}, fh, sort_keys=False)


if 'snakemake' in globals():
    record_config(json.loads(snakemake.params.context),
                  json.loads(snakemake.params.effective), snakemake.output.effective)
