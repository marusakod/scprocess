"""Exercise UMAP configuration/dispatch without a GPU or expression data."""
import argparse
import ast
import json
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, Mock, patch

import jsonschema
import numpy as np
import pytest

from scripts.scprocess_utils import _get_default_config_from_schema


ROOT = Path(__file__).resolve().parents[1]
TREE = ast.parse((ROOT / 'scripts/integration.py').read_text())


@pytest.fixture
def integration():
    # Execute the actual functions, omitting heavy module-level imports. All
    # expression-data IO and CPU/GPU backends are mocked below.
    namespace = {'argparse': argparse, 'np': np}
    functions = ast.Module(
        body=[n for n in TREE.body if isinstance(n, ast.FunctionDef)],
        type_ignores=[])
    exec(compile(functions, 'integration.py', 'exec'), namespace)
    return namespace


@pytest.mark.parametrize('workflow', ['config', 'zoom', 'join'])
def test_schema_defaults_and_validation(workflow):
    schema = json.loads((ROOT / f'resources/schemas/{workflow}.schema.json').read_text())
    jsonschema.Draft7Validator.check_schema(schema)
    for config in ({}, {'integration': {}}, {'integration': {'int_n_dims': 30}}):
        defaults = _get_default_config_from_schema(schema, config)
        assert defaults['integration']['int_umap_n_iter'] == 2000
    integration_schema = schema['properties']['integration']
    for value in (1, 750, 2000, 3500):
        config = {'integration': {'int_umap_n_iter': value}}
        effective = _get_default_config_from_schema(schema, config)['integration']
        effective.update(config['integration'])
        jsonschema.validate(effective, integration_schema)
        assert effective['int_umap_n_iter'] == value
    field = integration_schema['properties']['int_umap_n_iter']
    for value in (0, -1, 2.5, '2000', True, None):
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(value, field)


@pytest.mark.parametrize('command', ['run_integration', 'run_zoom_integration'])
@pytest.mark.parametrize('value', [None, '3500', '0', '-1', '2.5'])
def test_cli_default_override_and_rejection(integration, command, value):
    namespace = integration
    namespace.update(run_integration=Mock(), run_zoom_integration=Mock())
    main = next(n for n in TREE.body if isinstance(n, ast.If))
    body = ast.Module(body=main.body, type_ignores=[])
    argv = ['integration.py', command]
    if value is not None:
        argv += ['--umap_n_iter', value]
    scanpy = MagicMock()
    with patch.object(sys, 'argv', argv), patch.dict(sys.modules, {
            'scanpy': scanpy, 'scanpy.external': scanpy.external}):
        if value in ('0', '-1', '2.5'):
            with pytest.raises(SystemExit) as error:
                exec(compile(body, 'integration.py', 'exec'), namespace)
            assert error.value.code == 2
            namespace[command].assert_not_called()
        else:
            exec(compile(body, 'integration.py', 'exec'), namespace)
            assert namespace[command].call_args.kwargs['umap_n_iter'] == (
                2000 if value is None else int(value))


@pytest.mark.parametrize('gpu', [False, True])
@pytest.mark.parametrize('paga', [False, True])
@pytest.mark.parametrize('iterations', [None, 3500])
def test_backend_receives_iterations_and_preserves_initialization(integration, gpu, paga, iterations):
    adata = SimpleNamespace(obs={'sample_id': Mock()}, obsm={'X_pca': np.zeros((2, 2))})
    adata.obs['sample_id'].unique.return_value = ['s1', 's2']
    backend = MagicMock()
    integration.update(sc=backend, _get_clusts_from_adata=Mock(),
                       _get_embeddings_from_adata=Mock())
    kwargs = {} if iterations is None else {'umap_n_iter': iterations}
    scanpy = MagicMock()
    with patch.dict(sys.modules, {'scanpy': scanpy}):
        integration['_do_one_integration'](
            adata, 'sample_id', 'leiden', 2, [0.2], 'pca', gpu, 0.1,
            use_paga=paga, paga_cl='RNA_snn_res.0.2', **kwargs)
    backend.tl.umap.assert_called_once_with(
        adata, maxiter=2000 if iterations is None else iterations,
        init_pos='paga' if paga else 'spectral')


@pytest.mark.parametrize('workflow', ['standard', 'zoom', 'join_bpcells'])
def test_entry_points_forward_iterations(integration, workflow):
    adata = MagicMock(shape=(2, 2))
    pca_df = MagicMock(columns=['cell_id', 'pca_1', 'pca_2'])
    pca_df.select.return_value.to_numpy.return_value = np.zeros((2, 2))
    pca_df.__getitem__.return_value.to_list.return_value = ['c1', 'c2']
    integration.update(
        ad=SimpleNamespace(AnnData=Mock(return_value=adata)),
        pl=MagicMock(), gzip=MagicMock(), gc=Mock(),
        _get_hvg_mat=Mock(return_value=(MagicMock(), ['c1', 'c2'], [])),
        _get_cells_df=Mock(), _normalize_hvg_mat=Mock(),
        _calc_dbl_data=Mock(), _adata_filter_out_doublets=Mock(return_value=adata),
        _do_one_integration=Mock())
    integration['pl'].read_csv.return_value = pca_df
    common = dict(hvg_mat_f='counts.h5', sample_qc_f='qc.csv', coldata_f='cells.csv',
                  demux_type='none', exclude_mito=False, embedding='pca', n_dims=2,
                  cl_method='leiden', theta=0.1, res_ls_concat='0.2',
                  integration_f='unused.csv.gz', batch_var='sample_id',
                  umap_n_iter=3500)
    if workflow == 'standard':
        integration['run_integration'](
            **common, dbl_hvg_mat_f='doublets.h5', dbl_res=4, dbl_cl_prop=0.5)
    else:
        integration['run_zoom_integration'](
            **common, precomputed_pca_f='pca.csv.gz' if workflow == 'join_bpcells' else None)
    calls = integration['_do_one_integration'].call_args_list
    assert len(calls) == (2 if workflow == 'standard' else 1)
    assert all(call.kwargs['umap_n_iter'] == 3500 for call in calls)
    if workflow != 'standard':
        assert calls[0].kwargs['skip_pca'] == (workflow == 'join_bpcells')


@pytest.mark.parametrize('filename,rule', [
    ('integration.smk', 'run_integration'),
    ('zoom.smk', 'zoom_run_integration'),
    ('join.smk', 'join_integration')])
def test_rules_track_and_forward_config(filename, rule):
    text = (ROOT / 'rules' / filename).read_text()
    block = text.split(f'rule {rule}:', 1)[1].split('\nrule ', 1)[0]
    assert "['integration']['int_umap_n_iter']" in block
    assert '--umap_n_iter {params.umap_n_iter}' in block
    if filename == 'zoom.smk':
        assert "ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_umap_n_iter']" in block
