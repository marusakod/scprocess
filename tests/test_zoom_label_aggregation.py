import json
import pathlib
import sys

import jsonschema
import polars as pl
import pytest


SCRIPTS_DIR = pathlib.Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from label_celltypes import aggregate_predictions
from scprocess_utils import get_zoom_report_labels_f


SCHEMA_F = pathlib.Path(__file__).resolve().parents[1] / "resources/schemas/zoom.schema.json"


def test_reaggregate_preaggregated_labels_against_zoom_clusters(tmp_path):
    labels_f = tmp_path / "parent_labels.csv.gz"
    integration_f = tmp_path / "zoom_integration.csv.gz"

    pl.DataFrame({
        "model": ["test"] * 4,
        "sample_id": ["s1", "s1", "s2", "s2"],
        "cell_id": ["c1", "c2", "c3", "c4"],
        "hi_res_cl": ["parent1", "parent1", "parent2", "parent2"],
        "predicted_label_agg": ["old_A", "old_A", "old_B", "old_B"],
        "prop_hi_res_cl": [1.0, 1.0, 1.0, 1.0],
        "predicted_label_naive": ["A", "A", "B", "A"],
        "probability_naive": [0.9, 0.8, 0.7, 0.6],
    }).write_csv(labels_f)

    pl.DataFrame({
        "cell_id": ["c1", "c2", "c3", "c4"],
        "sample_id": ["s1", "s1", "s2", "s2"],
        "RNA_snn_res.2": ["zoom1", "zoom2", "zoom2", "zoom2"],
    }).write_csv(integration_f)

    result = aggregate_predictions(
        [str(labels_f)], str(integration_f), "RNA_snn_res.2", 0.5, "sample_id"
    ).sort("cell_id")

    assert result["hi_res_cl"].to_list() == ["zoom1", "zoom2", "zoom2", "zoom2"]
    assert result["predicted_label_agg"].to_list() == ["A", "A", "A", "A"]
    assert result["prop_hi_res_cl"].to_list() == [1.0, 2 / 3, 2 / 3, 2 / 3]


def test_save_cluster_names_requires_classifier_labels():
    with open(SCHEMA_F) as fh:
        schema = json.load(fh)
    zoom = {
        "zoom": {
            "name": "T_NK",
            "labels_source": "clusters",
            "labels_col": "RNA_snn_res.0.2",
            "sel_labels": ["0"],
            "save_cluster_names_file": True,
        }
    }

    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(zoom, schema)

    zoom["zoom"].update(labels_source="celltypist", model="Immune_All_Low")
    jsonschema.validate(zoom, schema)


def test_zoom_report_selects_aggregated_labels_for_classifier_sources():
    common = {"zoom": {"labels_source": "clusters"}}
    classifier = {"zoom": {"labels_source": "scprocess"}}

    assert get_zoom_report_labels_f(
        {"z": common}, "/project/zoom", "full", "2026-09-06", "z"
    ).endswith("/z/filtered_labels_full_z_2026-09-06.csv.gz")
    assert get_zoom_report_labels_f(
        {"z": classifier}, "/project/zoom", "full", "2026-09-06", "z"
    ).endswith("/z/aggregated_labels_full_z_2026-09-06.csv.gz")
