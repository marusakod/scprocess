import copy
import json
import os
from pathlib import Path
import subprocess
import sys

import polars as pl
import pytest
import yaml

from scripts.join_zoom import select_cells, write_csv
from scripts.scprocess_utils import get_join_zoom_specs, get_join_zoom_context, get_join_zoom_prepare_params

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def parent_fixture(tmp_path):
    root = tmp_path / 'joined'
    root.mkdir()
    scdata = tmp_path / 'scdata'
    scdata.mkdir()
    pl.DataFrame({'reference': ['mouse_2024'], 'gene_info_f': [str(scdata / 'genes.csv')]}).write_csv(scdata / 'index_parameters.csv')
    (scdata / 'genes.csv').touch()
    parent = {'join': {'name': 'atlas', 'proj_dir': str(root), 'date_stamp': '2026-09-15',
                      'ref_txome': 'mouse_2024', 'your_name': 'Test', 'affiliation': 'Test'},
              'projects': {}, 'zoom': ['zoom.yaml'],
              'integration': {'int_n_dims': 2, 'int_use_gpu': False, 'int_umap_n_iter': 2222},
              'hvg': {'hvg_n_hvgs': 10}, 'marker_genes': {'mkr_do_gsea': False},
              'shiny': {'app_tag': 'parent-app', 'annotation_csv': 'old-clusters.csv'}}
    manifest = {}
    for pid, sample in [('p_a', 's1'), ('p_b', 's2')]:
        source = tmp_path / pid
        source.mkdir()
        meta = source / 'metadata.csv'
        pl.DataFrame({'sample_id': [sample]}).write_csv(meta)
        cfg = {'project': {'proj_dir': str(source), 'full_tag': pid, 'short_tag': pid,
                          'date_stamp': '2026-09-15', 'sample_metadata': str(meta),
                          'fastq_dir': str(source / 'fastq'), 'ref_txome': 'mouse_2024',
                          'your_name': 'Test', 'affiliation': 'Test'}}
        cfg_f = source / 'config.yaml'
        cfg_f.write_text(yaml.safe_dump(cfg))
        parent['projects'][pid] = {'config': str(cfg_f)}
        source_int = source / 'output' / f'{pid}_integration'
        source_int.mkdir(parents=True)
        h5ad = source / f'{sample}.h5ad'
        h5ad.touch()
        (source_int / f'h5ads_clean_paths_{pid}_2026-09-15.yaml').write_text(yaml.safe_dump({sample: str(h5ad)}))
        manifest[f'{pid}_{sample}'] = {'path': str(h5ad), 'project_id': pid}
    out = root / 'output/atlas_join'
    out.mkdir(parents=True)
    cells = pl.DataFrame({'cell_id': [f'cell{i}' for i in range(8)],
                          'sample_id': ['p_a_s1']*4 + ['p_b_s2']*4,
                          'project_id': ['p_a']*4 + ['p_b']*4,
                          'RNA_snn_res.0.5': ['cl1']*3 + ['cl2'] + ['cl1']*3 + ['cl2'],
                          'UMAP1': list(range(8))})
    coldata = cells.with_columns(pl.lit('neurons').alias('cell_class'))
    meta = cells.select('sample_id', 'project_id').unique().with_columns(pl.lit('cortex').alias('tissue'))
    qc = cells.select('cell_id', 'sample_id', 'project_id').with_columns(
        pl.lit(3.0).alias('log_counts'), pl.lit(2.0).alias('log_feats'))
    for frame, stem, suffix in [(cells, 'integrated_dt', '.csv.gz'), (coldata, 'joint_coldata', '.csv.gz'),
                                 (meta, 'joint_sample_meta', '.csv'), (qc, 'joint_qc', '.csv.gz')]:
        write_csv(frame, out / f'{stem}_atlas_join_2026-09-15{suffix}')
    (out / 'h5ads_clean_paths_atlas_join_2026-09-15.yaml').write_text(yaml.safe_dump(manifest))
    spec = {'zoom': {'name': 'neurons', 'labels_source': 'clusters',
                     'labels_col': 'RNA_snn_res.0.5', 'sel_labels': ['cl1']},
            'qc': {'qc_min_cells': 2}}
    (root / 'zoom.yaml').write_text(yaml.safe_dump(spec))
    config_f = root / 'join.yaml'
    config_f.write_text(yaml.safe_dump(parent))
    return parent, spec, scdata, config_f, cells, coldata, meta, qc


def test_context_inherits_settings_and_isolates_outputs(parent_fixture):
    parent, spec, scdata, _, *_ = parent_fixture
    original = copy.deepcopy(parent)
    config, ctx = get_join_zoom_context(parent, 'neurons', ROOT, scdata)
    assert parent == original
    assert config['integration']['int_umap_n_iter'] == 2222
    assert config['integration']['int_use_gpu'] is False
    assert config['shiny']['app_tag'] == 'atlas_join_zoom_neurons'
    assert 'annotation_csv' not in config['shiny']
    assert ctx['out_dir'].endswith('/output/atlas_join_zoom/neurons')
    assert ctx['logs_dir'].endswith('/.log/join_zoom/neurons')
    assert not Path(ctx['out_dir']).exists()  # resolving/dry-running writes nothing
    assert all('/sources/' in p for p in ctx['project']['integrated_fs'])


@pytest.mark.parametrize('section,key,value,changes_preparation', [
    ('marker_genes', 'mkr_sel_res', 0.2, False),
    ('hvg', 'hvg_n_hvgs', 12, False),
    ('hvg', 'hvg_chunk_size', 100, False),
    ('resources', 'n_join_marker_genes_calc', 32, False),
    ('shiny', 'app_title', 'Updated title', False),
    ('integration', 'int_umap_n_iter', 3333, False),
    ('integration', 'int_n_dims', 3, True),
    ('integration', 'int_batch_var', 'project_id', True),
    ('qc', 'qc_min_cells', 3, True),
    ('zoom', 'sel_labels', ['cl2'], True),
])
def test_preparation_tracks_only_consumed_settings(parent_fixture, section, key, value, changes_preparation):
    parent, spec, scdata, config_f, *_ = parent_fixture
    cfg, ctx = get_join_zoom_context(parent, 'neurons', ROOT, scdata)
    before = get_join_zoom_prepare_params(cfg, ctx)
    spec.setdefault(section, {})[key] = value
    config_f.with_name('zoom.yaml').write_text(yaml.safe_dump(spec))
    cfg, ctx = get_join_zoom_context(parent, 'neurons', ROOT, scdata)
    after = get_join_zoom_prepare_params(cfg, ctx)
    assert (before != after) == changes_preparation


def test_preparation_resume_and_config_record_are_independent(parent_fixture, tmp_path):
    """Execute real preparation on eight cell headers, then edit the zoom YAML.

    No counts/PCA/marker calculations run. Input mtimes must still be honoured.
    """
    parent, spec, scdata, config_f, cells, *_ = parent_fixture
    _, ctx = get_join_zoom_context(parent, 'neurons', ROOT, scdata)
    hvg_python = os.environ.get('SCPROCESS_TEST_HVG_PYTHON', sys.executable)
    check = subprocess.run([hvg_python, '-c', 'import h5py, polars'], capture_output=True)
    if check.returncode:
        pytest.skip('Set SCPROCESS_TEST_HVG_PYTHON to an existing HVG environment')
    fixture_code = '''import h5py, json, sys
with h5py.File(sys.argv[1], 'w') as fh:
    fh.create_group('X')
    fh.create_group('obs').create_dataset('_index', data=json.loads(sys.argv[2]), dtype=h5py.string_dtype())
'''
    for batch in ctx['batches'].values():
        ids = cells.filter(pl.col('project_id') == batch['project_id'])['cell_id'].to_list()
        subprocess.run([hvg_python, '-c', fixture_code, batch['path'], json.dumps(ids)], check=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith('SNAKEMAKE_')}
    env.update(SCPROCESS_DATA_DIR=str(scdata), XDG_CACHE_HOME=str(tmp_path / 'cache'),
               PATH=str(Path(sys.executable).parent) + os.pathsep + env.get('PATH', ''),
               POLARS_MAX_THREADS='1', OPENBLAS_NUM_THREADS='1')
    wrapper = tmp_path / 'resume.smk'
    wrapper.write_text(f'include: {str(ROOT / "rules/join_zoom.smk")!r}\n'
        + (f'workflow.get_rule("join_zoom_prepare").conda_env = {str(Path(hvg_python).parent.parent)!r}\n'
           if hvg_python != sys.executable else ''))
    deployment = ['--software-deployment-method', 'conda'] if hvg_python != sys.executable else []
    command = [sys.executable, '-m', 'snakemake', '--snakefile', str(wrapper), *deployment,
        '--configfile', str(config_f), '--config', f'scprocess_dir={ROOT}', 'join_zoom_name=neurons',
        '--directory', str(config_f.parent), '--cores', '1', '--executor', 'local',
        '--profile', 'none', '--workflow-profile', 'none', '--rerun-triggers', 'mtime', 'params']

    def run(target, dry=False):
        args = command + (['--dry-run'] if dry else []) + ['--', target]
        proc = subprocess.run(args, capture_output=True, text=True, timeout=45, env=env)
        output = proc.stdout + proc.stderr
        assert proc.returncode == 0, output
        return output

    run('join_zoom_record_config')
    assert 'Nothing to be done' in run('join_zoom_record_config', dry=True)
    prepared = [Path(ctx['selected_f']), *map(Path, ctx['project']['integrated_fs'])]
    before = {str(p): (p.stat().st_mtime_ns, p.read_bytes()) for p in prepared}
    # Comments/YAML formatting alone must not invalidate either rule.
    zoom_f = config_f.with_name('zoom.yaml')
    zoom_f.write_text(zoom_f.read_text() + '\n# No analytical change\n')
    assert 'Nothing to be done' in run('join_zoom_record_config', dry=True)
    spec['marker_genes'] = {'mkr_sel_res': 0.2}
    spec['resources'] = {'n_join_marker_genes_calc': 32}
    zoom_f.write_text(yaml.safe_dump(spec))
    assert 'Nothing to be done' in run('join_zoom_prepare', dry=True)
    output = run('join_zoom_record_config', dry=True)
    assert 'rule join_zoom_record_config:' in output
    assert 'rule join_zoom_prepare:' not in output
    run('join_zoom_record_config')
    assert before == {str(p): (p.stat().st_mtime_ns, p.read_bytes()) for p in prepared}
    recorded = yaml.safe_load(Path(ctx['effective_f']).read_text())
    assert recorded['config']['marker_genes']['mkr_sel_res'] == 0.2
    spec['qc']['qc_min_cells'] = 3
    zoom_f.write_text(yaml.safe_dump(spec))
    assert 'rule join_zoom_prepare:' in run('join_zoom_prepare', dry=True)
    spec['qc']['qc_min_cells'] = 2
    zoom_f.write_text(yaml.safe_dump(spec))
    assert 'Nothing to be done' in run('join_zoom_prepare', dry=True)
    # Genuine input updates are never hidden using ancient() or ignored mtimes.
    h5ad = Path(next(iter(ctx['batches'].values()))['path'])
    # Snakemake also checks contents for small inputs: change this tiny fixture,
    # not just its timestamp, to represent a genuinely updated source matrix.
    subprocess.run([hvg_python, '-c',
        'import h5py,sys\nwith h5py.File(sys.argv[1], "r+") as fh: fh.attrs["revision"] = 2',
        str(h5ad)], check=True)
    newer = max(p.stat().st_mtime_ns for p in prepared) + 2_000_000_000
    os.utime(h5ad, ns=(newer, newer))
    assert 'rule join_zoom_prepare:' in run('join_zoom_prepare', dry=True)


@pytest.mark.parametrize('change', [
    {'cell_cycle': {'cyc_target_cells': 100}}, {'train_xgboost': {}},
    {'integration': {'int_cell_cycle_correction': 'within_sample'}},
    {'zoom': {}}, {'hvg': {'hvg_method': 'all'}},
])
def test_unsupported_or_empty_settings_fail(parent_fixture, change):
    parent, spec, _, config_f, *_ = parent_fixture
    spec.update(change)
    config_f.with_name('zoom.yaml').write_text(yaml.safe_dump(spec))
    with pytest.raises(ValueError, match='Invalid joined zoom'):
        get_join_zoom_specs(parent, ROOT)


def test_duplicate_names_and_bad_selector(parent_fixture):
    parent, _, scdata, config_f, *_ = parent_fixture
    with pytest.raises(ValueError, match='Unknown joined zoom'):
        get_join_zoom_context(parent, 'missing', ROOT, scdata)
    parent['zoom'].append('zoom.yaml')
    with pytest.raises(ValueError, match='Duplicate'):
        get_join_zoom_specs(parent, ROOT)


def test_disabled_paga_and_partial_overrides(parent_fixture):
    parent, spec, scdata, config_f, *_ = parent_fixture
    spec['integration'] = {'int_use_paga': False, 'int_res_ls': [0.5], 'int_umap_n_iter': 2000}
    spec['marker_genes'] = {'mkr_sel_res': 0.5}
    config_f.with_name('zoom.yaml').write_text(yaml.safe_dump(spec))
    config, _ = get_join_zoom_context(parent, 'neurons', ROOT, scdata)
    assert not config['integration']['int_use_paga']
    assert 'int_paga_cl_res' not in config['integration']
    assert config['integration']['int_umap_n_iter'] == 2000
    assert config['integration']['int_n_dims'] == 2


def test_mismatched_paga_resolution_is_rejected(parent_fixture):
    parent, spec, scdata, config_f, *_ = parent_fixture
    spec['integration'] = {'int_res_ls': [0.2, 0.5]}
    config_f.with_name('zoom.yaml').write_text(yaml.safe_dump(spec))
    with pytest.raises(ValueError, match='int_paga_cl_res'):
        get_join_zoom_context(parent, 'neurons', ROOT, scdata)


def test_app_tag_collision_is_rejected(parent_fixture):
    parent, spec, _, config_f, *_ = parent_fixture
    spec['shiny'] = {'app_tag': 'parent-app'}
    config_f.with_name('zoom.yaml').write_text(yaml.safe_dump(spec))
    with pytest.raises(ValueError, match='shares shiny.app_tag'):
        get_join_zoom_specs(parent, ROOT)


def test_selection_preserves_metadata_and_parent_coordinates(parent_fixture):
    _, spec, _, _, cells, coldata, meta, qc = parent_fixture
    result, summary = select_cells(cells, coldata, meta, cells, qc, spec['zoom'], spec['qc'], 2)
    assert result['cell_id'].to_list() == ['cell0', 'cell1', 'cell2', 'cell4', 'cell5', 'cell6']
    assert result['sample_id'].unique().sort().to_list() == ['p_a_s1', 'p_b_s2']
    assert result['tissue'].unique().to_list() == ['cortex']
    assert 'parent_RNA_snn_res.0.5' in result.columns and 'parent_UMAP1' in result.columns
    assert summary['n_retained'].sum() == 6


def test_selection_can_retain_one_project(parent_fixture):
    _, spec, _, _, cells, coldata, meta, qc = parent_fixture
    spec['zoom']['exclude'] = {'sample_id': ['p_b_s2']}
    result, summary = select_cells(cells, coldata, meta, cells, qc, spec['zoom'], spec['qc'], 2)
    assert result['project_id'].unique().to_list() == ['p_a']
    assert summary['n_retained'].to_list() == [3, 0]


@pytest.mark.parametrize('failure', ['duplicate', 'missing_label', 'empty', 'missing_qc'])
def test_selection_errors(parent_fixture, failure):
    _, spec, _, _, cells, coldata, meta, qc = parent_fixture
    labels = cells
    if failure == 'duplicate':
        labels = pl.concat([cells, cells.head(1)])
    elif failure == 'missing_label':
        spec['zoom']['sel_labels'] = ['absent']
    elif failure == 'empty':
        spec['qc']['qc_min_counts'] = 100000
    else:
        qc = qc.tail(7)
    with pytest.raises(ValueError):
        select_cells(cells, coldata, meta, labels, qc, spec['zoom'], spec['qc'], 2)


@pytest.mark.parametrize('gpu,backend', [(False, 'bpcells'), (True, 'bpcells'), (False, 'scanpy')])
def test_real_join_zoom_dag(parent_fixture, tmp_path, gpu, backend):
    parent, spec, scdata, config_f, *_ = parent_fixture
    parent['integration'].update(int_use_gpu=gpu, int_pca_method=backend)
    config_f.write_text(yaml.safe_dump(parent))
    env = {k: v for k, v in os.environ.items() if not k.startswith('SNAKEMAKE_')}
    env.update(SCPROCESS_DATA_DIR=str(scdata), XDG_CACHE_HOME=str(tmp_path / 'cache'))
    profile = tmp_path / 'profile'
    profile.mkdir()
    (profile / 'config.yaml').write_text(yaml.safe_dump({'set-resources': {
        'join_integration': {'gpu': 1 if gpu else 0, 'slurm_partition': 'test_gpu' if gpu else 'test_cpu'}}}))
    proc = subprocess.run([sys.executable, '-m', 'snakemake', '--snakefile', str(ROOT / 'rules/join_zoom.smk'),
        '--configfile', str(config_f), '--config', f'scprocess_dir={ROOT}', 'join_zoom_name=neurons',
        '--directory', str(config_f.parent), '--dry-run', '--printshellcmds', '--cores', '2',
        '--executor', 'local', '--profile', 'none', '--workflow-profile', str(profile), '--nolock', 'all'],
        capture_output=True, text=True, timeout=45, env=env)
    output = proc.stdout + proc.stderr
    assert proc.returncode == 0, output
    for rule in ('join_zoom_prepare', 'join_zoom_hvg_stats', 'join_zoom_hvg_trend',
                 'join_zoom_hvg_clip', 'join_zoom_hvgs', 'join_integration', 'join_make_pseudobulks',
                 'join_marker_genes_prep', 'join_marker_genes_calc', 'join_render_html'):
        assert f'rule {rule}:' in output, output
    assert '--retain-metadata' in output and '--require-all-cells' in output
    assert '--umap_n_iter 2222' in output
    integration_block = output.split('rule join_integration:', 1)[1].split('Shell command:', 1)[0]
    assert f'gpu={int(gpu)}' in integration_block
    # String resources are evaluated later when the generated inputs exist.
    assert 'slurm_partition=' in integration_block
    assert ('rule join_pca:' in output) == (backend == 'bpcells')
    for rule in ('run_integration', 'make_clean_h5ads', 'run_cellbender', 'create_doublets_hvg_matrix'):
        assert f'rule {rule}:' not in output
    assert not (config_f.parent / 'output/atlas_join_zoom').exists()


@pytest.mark.parametrize('zoom', [False, True])
def test_zoom_only_params_do_not_invalidate_standard_joins(parent_fixture, tmp_path, zoom):
    """Check actual parsed rule parameters, including absence of no-op values."""
    parent, _, scdata, config_f, *_ = parent_fixture
    parent['shiny'].pop('annotation_csv')
    config_f.write_text(yaml.safe_dump(parent))
    workflow = ROOT / 'rules' / ('join_zoom.smk' if zoom else 'join.smk')
    wrapper = tmp_path / 'inspect_params.smk'
    wrapper.write_text(
        f'include: {str(workflow)!r}\n'
        'print("JOIN_PARAMS=" + json.dumps({name: dict(workflow.get_rule(name).params.items()) '
        'for name in ("join_build_matrix", "join_integration")}))\n'
    )
    env = {k: v for k, v in os.environ.items() if not k.startswith('SNAKEMAKE_')}
    env.update(SCPROCESS_DATA_DIR=str(scdata), XDG_CACHE_HOME=str(tmp_path / 'cache'))
    config_args = [f'scprocess_dir={ROOT}']
    if zoom:
        config_args.append('join_zoom_name=neurons')
    proc = subprocess.run([sys.executable, '-m', 'snakemake', '--snakefile', str(wrapper),
        '--configfile', str(config_f), '--config', *config_args,
        '--directory', str(config_f.parent), '--list-rules', '--cores', '1',
        '--executor', 'local', '--profile', 'none', '--workflow-profile', 'none'],
        capture_output=True, text=True, timeout=45, env=env)
    output = proc.stdout + proc.stderr
    assert proc.returncode == 0, output
    params = json.loads(next(line.removeprefix('JOIN_PARAMS=') for line in output.splitlines()
                             if line.startswith('JOIN_PARAMS=')))
    matrix = params['join_build_matrix']
    integration = params['join_integration']
    if zoom:
        assert matrix.pop('require_all_cells') == '--require-all-cells'
        assert integration.pop('retain_metadata') == '--retain-metadata'
    assert matrix == {'project_ids': 'p_a p_b'}
    assert set(integration) == {'embedding', 'n_dims', 'umap_n_iter', 'cl_method',
        'theta_concat', 'batch_var_concat', 'res_ls_concat', 'use_paga', 'paga_cl_res',
        'int_use_gpu', 'pca_method'}


@pytest.mark.parametrize('zoom', [False, True])
def test_h5ad_manifest_rule_tracks_retained_batches(parent_fixture, tmp_path, zoom):
    """Execute the real rule, including Snakemake's output validation/resume."""
    parent, _, scdata, config_f, *_ = parent_fixture
    parent['shiny'].pop('annotation_csv')
    config_f.write_text(yaml.safe_dump(parent))
    parent_manifest_f = config_f.parent / 'output/atlas_join/h5ads_clean_paths_atlas_join_2026-09-15.yaml'
    original_manifest = parent_manifest_f.read_bytes()
    sources = yaml.safe_load(original_manifest)
    if zoom:
        _, ctx = get_join_zoom_context(parent, 'neurons', ROOT, scdata)
        out_dir, tag = Path(ctx['out_dir']), ctx['tag']
        manifests = [Path(p) for p in ctx['project']['h5ads_yaml_fs']]
        for i, path in enumerate(manifests):
            path.parent.mkdir(parents=True, exist_ok=True)
            # All cells from p_b were excluded by the zoom selection/QC.
            path.write_text(yaml.safe_dump({'s1': sources['p_a_s1']} if i == 0 else {}))
        expected = {'p_a_s1'}
    else:
        out_dir, tag = config_f.parent / 'output/atlas_join', 'atlas_join'
        expected = {'p_a_s1', 'p_b_s2'}
    output_manifest = out_dir / f'h5ads_clean_paths_{tag}_2026-09-15.yaml'
    if not zoom:
        # The fixture's parent manifest is an input only for zooms. A standard
        # join should create this output from its two source-project manifests.
        output_manifest.unlink()
    env = {k: v for k, v in os.environ.items() if not k.startswith('SNAKEMAKE_')}
    env.update(SCPROCESS_DATA_DIR=str(scdata), XDG_CACHE_HOME=str(tmp_path / 'cache'),
               PATH=str(Path(sys.executable).parent) + os.pathsep + env.get('PATH', ''))
    config_args = [f'scprocess_dir={ROOT}']
    if zoom:
        config_args.append('join_zoom_name=neurons')
    command = [sys.executable, '-m', 'snakemake', '--snakefile',
        str(ROOT / 'rules' / ('join_zoom.smk' if zoom else 'join.smk')),
        '--configfile', str(config_f), '--config', *config_args,
        '--directory', str(config_f.parent), '--cores', '1', '--executor', 'local',
        '--profile', 'none', '--workflow-profile', 'none', '--rerun-triggers', 'mtime', 'params',
        '--allowed-rules', 'join_build_h5ads_yaml', '--', 'join_build_h5ads_yaml']

    def run(extra=()):
        # Insert options before the end-of-options marker.
        args = command[:-2] + list(extra) + command[-2:]
        proc = subprocess.run(args, capture_output=True, text=True, timeout=45, env=env)
        output = proc.stdout + proc.stderr
        assert proc.returncode == 0, output
        return output

    run()
    actual = yaml.safe_load(output_manifest.read_text())
    assert set(actual) == expected
    for key, entry in actual.items():
        link = Path(entry['path'])
        assert link.is_symlink() and link.resolve() == Path(sources[key]['path']).resolve()
    if zoom:
        assert (out_dir / 'h5ads/.snakemake_timestamp').exists()
        assert not (out_dir / 'h5ads/p_b_s2.h5ad').exists()
    else:
        assert not (out_dir / 'h5ads/.snakemake_timestamp').exists()
    assert 'Nothing to be done' in run(['--dry-run'])
    if zoom:
        # A changed subset must remove stale links, while preserving sources.
        manifests[0].write_text('{}\n')
        manifests[1].write_text(yaml.safe_dump({'s2': sources['p_b_s2']}))
        run(['--forcerun', 'join_build_h5ads_yaml'])
        assert set(yaml.safe_load(output_manifest.read_text())) == {'p_b_s2'}
        assert not (out_dir / 'h5ads/p_a_s1.h5ad').exists()
        assert (out_dir / 'h5ads/p_b_s2.h5ad').is_file()
        assert parent_manifest_f.read_bytes() == original_manifest
    assert all(Path(entry['path']).is_file() for entry in sources.values())


def test_join_zoom_shiny_paths(parent_fixture, tmp_path):
    parent, _, scdata, config_f, *_ = parent_fixture
    cfg, ctx = get_join_zoom_context(parent, 'neurons', ROOT, scdata)
    out = Path(ctx['out_dir'])
    out.mkdir(parents=True)
    tag = ctx['tag']
    for stem, suffix in [('h5ads_clean_paths', '.yaml'), ('joint_sample_meta', '.csv'),
                         ('integrated_dt', '.csv.gz'), ('pb_marker_genes', '.csv.gz'), ('pb_hvgs', '.csv.gz')]:
        res = '_0.2' if stem.startswith('pb_') else ''
        (out / f'{stem}_{tag}{res}_2026-09-15{suffix}').touch()
    env = {k: v for k, v in os.environ.items() if not k.startswith('SNAKEMAKE_')}
    env.update(SCPROCESS_DATA_DIR=str(scdata), XDG_CACHE_HOME=str(tmp_path / 'cache'))
    proc = subprocess.run([sys.executable, '-m', 'snakemake', '--snakefile', str(ROOT / 'rules/shiny.smk'),
        '--configfile', str(config_f), '--config', f'scprocess_dir={ROOT}', 'join_zoom_name=neurons',
        '--directory', str(config_f.parent), '--dry-run', '--printshellcmds', '--cores', '2',
        '--executor', 'local', '--profile', 'none', '--workflow-profile', 'none', '--nolock', 'configure_shiny_app'],
        capture_output=True, text=True, timeout=45, env=env)
    output = proc.stdout + proc.stderr
    assert proc.returncode == 0, output
    assert 'public/shiny_atlas_join_zoom_neurons/' in output
    assert 'public/shiny_parent-app/' not in output
