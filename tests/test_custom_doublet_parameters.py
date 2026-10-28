import json
from pathlib import Path
import pytest
import jsonschema
from scripts import scprocess_utils as utils
from tests.test_qc_hvg_workflow_dependencies import setup_rules, run_snakemake

ROOT = Path(__file__).resolve().parents[1]

@pytest.mark.parametrize("level", ["sample_id", "pool_id"])
@pytest.mark.parametrize("value", [None, 0, 1])
def test_custom_schema_accepts_doublet_uncertainty(level, value):
    schema = json.loads((ROOT / "resources/schemas/custom_sample_params.schema.json").read_text())
    jsonschema.validate({level: {"run1": {"qc": {"dbl_dbr_sd": value}}}}, schema)

@pytest.mark.parametrize("value", [-1, "1", True])
def test_custom_schema_rejects_invalid_uncertainty(value):
    schema = json.loads((ROOT / "resources/schemas/custom_sample_params.schema.json").read_text())
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate({"sample_id": {"s1": {"qc": {"dbl_dbr_sd": value}}}}, schema)

@pytest.mark.parametrize("global_value,custom,expected", [
    (None, {}, None), (1, {}, 1),
    (None, {"s1": {"qc": {"dbl_dbr_sd": 1}}}, 1),
    (1, {"s1": {"qc": {"dbl_dbr_sd": None}}}, None),
    (1, {"s1": {"qc": {"dbl_dbr_sd": 0}}}, 0),
    (None, {"s2": {"qc": {"dbl_dbr_sd": 1}}}, None),
])
def test_run_resolver_preserves_override_and_inheritance(global_value, custom, expected):
    config = {"qc": {"dbl_dbr_sd": global_value}, "ambient": {"ambient_method": "none"}}
    result = utils._get_run_parameters_one_run("s1", config, None, custom)
    assert result["qc"]["dbl_dbr_sd"] == expected


def test_pooled_doublet_override_rejects_sample_level():
    from unittest.mock import patch
    config = {"multiplexing": {"demux_type": "hto"}}
    with patch.object(utils, "_get_custom_parameters", return_value={"s1": {"qc": {"dbl_dbr_sd": 1}}}):
        with pytest.raises(ValueError, match="must be configured under pool_id"):
            utils.get_run_parameters(config, None, "pool_id", ["p1"])

@pytest.mark.parametrize("value, rendered", [(1, "1"), (None, "NULL"), (0, "0")])
def test_qc_rule_passes_resolved_run_value(tmp_path, value, rendered):
    setup_rules(tmp_path, "qc.smk")
    snake = tmp_path / "Snakefile"
    text = snake.read_text().replace("config['qc']['dbl_dbr_sd']", repr(value))
    snake.write_text(text)
    for name in ["af/af_s1/af_counts_mat.h5", "qc/qc_thresholds_by_sample_id_test_date.csv",
                 "ambient/ambient_run_statistics_test_date.csv"]:
        p = tmp_path / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.touch()
    output = run_snakemake(tmp_path, "doublets/dbl_s1/scDblFinder_s1_outputs_test_date.csv.gz",
                           "--dry-run", "--printshellcmds")
    assert f"dbl_dbr_sd       = {rendered}" in output


def test_pooled_run_accepts_mixed_pool_and_sample_qc(tmp_path):
    custom = tmp_path / "custom.yaml"
    custom.write_text("pool_id:\n  p1:\n    qc:\n      dbl_dbr_sd: 1\nsample_id:\n  s1:\n    qc:\n      qc_min_counts: 800\n")
    config = {"project": {"is_flex": False, "custom_sample_params": str(custom)},
              "multiplexing": {"demux_type": "hto"},
              "ambient": {"ambient_method": "none"}, "qc": {"dbl_dbr_sd": None}}
    params, run_var = utils.get_run_parameters(config, None, "pool_id", ["p1"])
    assert run_var == "pool_id"
    assert params["p1"]["qc"]["dbl_dbr_sd"] == 1


def test_default_resolution_preserves_method_validation():
    import json
    from pathlib import Path
    import jsonschema
    from scripts.scprocess_utils import _get_default_config_from_schema
    from snakemake.utils import update_config
    schema = json.loads((Path(__file__).resolve().parents[1] / 'resources/schemas/config.schema.json').read_text())
    for method in (None, 'adaptive', 'manual'):
        supplied = {'integration': {}} if method is None else {'integration': {'int_dbl_threshold_method': method}}
        effective = _get_default_config_from_schema(schema, supplied)['integration']
        assert 'int_dbl_cl_prop' not in effective
        if method == 'manual':
            supplied['integration']['int_dbl_cl_prop'] = .3
        update_config(effective, supplied['integration'])
        jsonschema.validate(effective, schema['properties']['integration'])
    supplied = {'integration': {'int_dbl_threshold_method': 'adaptive', 'int_dbl_cl_prop': .3}}
    effective = _get_default_config_from_schema(schema, supplied)['integration']
    update_config(effective, supplied['integration'])
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(effective, schema['properties']['integration'])
