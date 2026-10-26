"""Excluded batch placeholders must not prevent full-cohort evaluation."""
from pathlib import Path
import sys
import pandas as pd
import polars as pl
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from label_celltypes import run_xgboost
from train_xgboost_workflow import aggregate_plan_predictions


def test_empty_placeholder_needs_no_model(tmp_path):
    path = tmp_path / 'excluded.h5ad'
    path.touch()
    result = run_xgboost('excluded', 'sample_id', 'model', path,
        'absent-model', 'absent-classes', 'absent-genes')
    assert result.height == 0
    assert result.schema['probability_naive'] == pl.Float64
    assert set(result.columns) == {'cell_id', 'predicted_label_naive',
        'probability_naive', 'labeller', 'sample_id', 'model', 'model_id', 'model_sha256'}
    with pytest.raises(FileNotFoundError):
        run_xgboost('missing', 'sample_id', 'model', tmp_path / 'missing.h5ad',
            'absent-model', 'absent-classes', 'absent-genes')


def test_header_only_batch_and_complete_membership(tmp_path):
    plan = tmp_path / 'plan.csv.gz'
    pd.DataFrame(dict(cell_id=['a', 'b'], label=['x', 'y'],
        original_label=['x', 'y'], split=['train', 'test'],
        selected=[True, False], label_supported_in_training=[True, True])).to_csv(plan, index=False)
    pred = tmp_path / 'pred.csv.gz'
    pd.DataFrame(dict(cell_id=['a', 'b'], predicted_label_naive=['x', 'x'],
        probability_naive=[.9, .8])).to_csv(pred, index=False)
    empty = tmp_path / 'empty.csv.gz'
    pd.DataFrame(columns=['cell_id', 'predicted_label_naive', 'probability_naive']).to_csv(empty, index=False)
    integration = tmp_path / 'integration.csv'
    pd.DataFrame({'cell_id': ['a', 'b']}).to_csv(integration, index=False)
    output = tmp_path / 'demo_fulldata_predictions.csv.gz'
    aggregate_plan_predictions([empty, pred], plan, output, integration)
    result = pd.read_csv(output)
    assert set(result.cell_id) == {'a', 'b'}
    assert result.probability.tolist() == [.9, .8]
    metrics = pd.read_csv(tmp_path / 'demo_evaluation.csv').set_index('split')
    assert metrics.loc['test', 'accuracy'] == 0
    pd.DataFrame({'cell_id': ['a', 'b', 'missing']}).to_csv(integration, index=False)
    with pytest.raises(ValueError, match='Integrated cells missing'):
        aggregate_plan_predictions([empty, pred], plan, output, integration)
