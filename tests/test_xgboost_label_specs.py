"""Label mapping, display specifications and local installation contracts."""
from pathlib import Path
import csv
import importlib.util
import json
import sys
import pandas as pd
import pytest
import yaml

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from train_xgboost import read_label_specs, label_mapping, get_xgboost_model_info
from train_xgboost import write_label_specifications
from scprocess_utils import _check_xgboost_label_files


def tables(tmp_path):
    fine=tmp_path/'labels.csv'
    fine.write_text('label,colour,coarse_label\nb,#223344,group2\na,#112233,group1\n')
    coarse=tmp_path/'coarse.csv'
    coarse.write_text('coarse_label,coarse_colour\ngroup1,#123456\ngroup2,#654321\n')
    return fine,coarse


def setup_module_api():
    spec=importlib.util.spec_from_file_location('scprocess_setup_specs_test',ROOT/'scripts/setup.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def test_default_has_no_coarse_labels():
    fine,coarse=read_label_specs(labels=['b','a'])
    assert [row['label'] for row in fine]==['b','a']
    assert coarse==[] and label_mapping(fine)=={}


def test_independent_order_and_colours(tmp_path):
    fine,coarse=tables(tmp_path)
    labels,groups=read_label_specs(fine,coarse,['a','b'])
    assert [row['label'] for row in labels]==['b','a']
    assert [row['coarse_label'] for row in groups]==['group1','group2']
    assert labels[0]['colour']=='#223344'
    assert label_mapping(labels)=={'b':'group2','a':'group1'}


@pytest.mark.parametrize('alias',['annotation','fine_label','label'])
def test_legacy_input_names(tmp_path,alias):
    path=tmp_path/'labels.csv';path.write_text(f'{alias},coarse_label\na,A\n')
    labels,_=read_label_specs(label_map_f=path,labels=['a'])
    assert labels[0]['label']=='a'


@pytest.mark.parametrize('contents,match',[
    ('label\na\na\n','duplicate'),
    ('label,colour\na,red\n','Invalid colour'),
    ('label,colour\na,#112233\nb,\n','every label'),
    ('label,coarse_label\na,A\nb,\n','every label'),
    ('label,annotation\na,a\n','one label column')])
def test_invalid_specifications_fail(tmp_path,contents,match):
    path=tmp_path/'labels.csv';path.write_text(contents)
    with pytest.raises(ValueError,match=match):read_label_specs(path)


def test_mapping_coverage_and_coarse_groups(tmp_path):
    fine,coarse=tables(tmp_path)
    with pytest.raises(ValueError,match='miss model labels'):read_label_specs(fine,coarse,['unknown'])
    coarse.write_text('coarse_label\ngroup1\n')
    with pytest.raises(ValueError,match='miss groups'):read_label_specs(fine,coarse)
    with pytest.raises(ValueError,match='not both'):read_label_specs(fine,label_map_f=fine)


def test_display_only_and_optional_coarse_table(tmp_path):
    fine,coarse=tables(tmp_path)
    labels,groups=read_label_specs(fine)
    assert [row['coarse_label'] for row in groups]==['group2','group1']
    fine.write_text('label,colour\nb,#223344\na,#112233\n')
    assert read_label_specs(fine)[1]==[]
    with pytest.raises(ValueError,match='require a coarse_label mapping'):read_label_specs(fine,coarse)


def test_relative_paths_remain_project_relative(tmp_path):
    directory=tmp_path/'data/configs';directory.mkdir(parents=True)
    (directory/'cells.csv').write_text('cell_id,label\nc,a\n')
    config={'annots_f':'data/configs/cells.csv'}
    _check_xgboost_label_files(config,tmp_path)
    assert config['annots_f']==str(directory/'cells.csv')
    # A bare filename resolves at the project root, not data/configs.
    with pytest.raises(FileNotFoundError):_check_xgboost_label_files({'annots_f':'cells.csv'},tmp_path)


def test_export_and_coarse_metrics_use_direct_mapping(tmp_path):
    fine=tmp_path/'fine.csv';fine.write_text('label,coarse_label\nb,Neuron\na,Neuron\n')
    classes=tmp_path/'classes.csv';classes.write_text('class\na\nb\n')
    counts=tmp_path/'counts.csv';counts.write_text('label,len\na,2\nb,2\n')
    predictions=tmp_path/'predictions.csv.gz'
    pd.DataFrame(dict(label=['a','b'],predicted_label=['b','a'],split=['test','test'],
        selected=[False,False])).to_csv(predictions,index=False)
    write_label_specifications(classes, counts, predictions, tmp_path, xgboost_tag='demo', label_specs_f=str(fine), coarse_label_specs_f=None, label_map_f=None)
    result=pd.read_csv(tmp_path/'demo_coarse_evaluation.csv')
    assert result.loc[0,'accuracy']==1 and result.loc[0,'n_cells']==2
    assert pd.read_csv(tmp_path/'demo_label_specs.csv').label.tolist()==['b','a']
    assert set(pd.read_csv(predictions).predicted_label)=={'a','b'}


def test_no_coarse_mapping_writes_empty_metrics(tmp_path):
    classes=tmp_path/'classes.csv';classes.write_text('class\na\n')
    counts=tmp_path/'counts.csv';counts.write_text('label,len\na,1\n')
    write_label_specifications(classes, counts, 'unused', tmp_path, xgboost_tag='demo', label_specs_f=None, coarse_label_specs_f=None, label_map_f=None)
    assert pd.read_csv(tmp_path/'demo_coarse_evaluation.csv').empty
    assert pd.read_csv(tmp_path/'demo_coarse_label_specs.csv').empty


def model_source(tmp_path):
    source=tmp_path/'source';source.mkdir()
    (source/'tag_xgboost_model.json').write_text('weights')
    (source/'tag_allowed_cls.csv').write_text('class\na\nb\n')
    (source/'tag_selected_genes.txt').write_text('g1\n')
    fine,coarse=tables(tmp_path)
    (source/'tag_label_specs.csv').write_bytes(fine.read_bytes())
    (source/'tag_coarse_label_specs.csv').write_bytes(coarse.read_bytes())
    return source


def test_installation_preserves_specifications_and_checksums(tmp_path):
    setup=setup_module_api();source=model_source(tmp_path);data=tmp_path/'data'
    setup.add_xgboost_classifier(str(data),'local','tag',str(source),False)
    dest=data/'xgboost/local'
    assert (dest/'local_label_specs.csv').read_bytes()==(source/'tag_label_specs.csv').read_bytes()
    assert (dest/'local_coarse_label_specs.csv').read_bytes()==(source/'tag_coarse_label_specs.csv').read_bytes()
    setup.save_available_xgboost_csv(str(data/'xgboost/xgboost_models.csv'),str(data),['local'])
    row=pd.read_csv(data/'xgboost/xgboost_models.csv').iloc[0]
    assert row.label_specs_f==str(dest/'local_label_specs.csv')
    get_xgboost_model_info(dest/'local_xgboost_model.json','local')
    (dest/'local_label_specs.csv').write_text('changed')
    with pytest.raises(ValueError,match='checksum'):get_xgboost_model_info(dest/'local_xgboost_model.json','local')


def test_setup_direct_source_and_current_config_tag(tmp_path):
    setup=setup_module_api();source=model_source(tmp_path)
    assert setup._resolve_xgboost_spec({'src_dir':str(source),'name':'custom'})==('custom','tag',str(source))
    config=tmp_path/'config.yaml'
    config.write_text(yaml.safe_dump({'project':{'proj_dir':str(tmp_path),'short_tag':'short','full_tag':'full'},
        'train_xgboost':{'xgboost_tag':'new'}}))
    name,tag,path=setup._resolve_xgboost_spec({'config':str(config),'name':'custom'})
    assert (name,tag,path)==('custom','xgboost_new',str(tmp_path/'output/short_train_xgboost'))


def test_schemas_accept_optional_specifications():
    from jsonschema import Draft202012Validator
    for name in ['config','join','zoom']:
        schema=json.loads((ROOT/f'resources/schemas/{name}.schema.json').read_text())['properties']['train_xgboost']
        Draft202012Validator(schema).validate({'annots_f':'data/configs/cells.csv',
            'label_specs_f':'data/configs/labels.csv','coarse_label_specs_f':'data/configs/coarse.csv'})



def test_cli_optional_features_can_be_omitted(tmp_path):
  import subprocess
  import polars as pl

  rows = [dict(cell_id = f'{batch}_{label}_{cell}', label = label, sample_id = batch)
    for batch in ['s1', 's2', 's3', 's4'] for label in ['a', 'b'] for cell in range(5)]
  data = pl.DataFrame(rows)
  annotations = tmp_path / 'annotations.csv'
  integration = tmp_path / 'integration.csv'
  manifest = tmp_path / 'h5ads.yaml'
  data.select('cell_id', 'label').write_csv(annotations)
  data.select('cell_id', 'sample_id').write_csv(integration)
  manifest.write_text(yaml.safe_dump({batch: 'unused.h5ad' for batch in data['sample_id'].unique()}))
  training_metadata_f = tmp_path / 'training_metadata.csv.gz'
  settings = dict(annots_f = annotations, integration_f = integration, h5ads_yaml = manifest,
    batch_var = 'sample_id', training_metadata_f = training_metadata_f, coverage_f = tmp_path / 'coverage.csv',
    counts_f = tmp_path / 'counts.csv', validation_group_var = 'batch', is_join = 0,
    train_cells_per_label = 20, validation_cells_per_label = 10, min_cells_per_label = 2,
    seed = 7, cluster_majority_vote_min_purity = .9, validation_fraction = .25)
  command = [sys.executable, str(ROOT / 'scripts/train_xgboost.py'), 'prepare']
  for name, value in settings.items():
    command += ['--' + name, str(value)]
  subprocess.run(command, check = True, capture_output = True, text = True, timeout = 30)
  plan = pl.read_csv(training_metadata_f)
  assert set(plan['split']) == {'train', 'validation'}
  assert plan['label'].equals(plan['original_label'])

  # Coupled test settings must still be supplied together.
  invalid = subprocess.run(command + ['--test_group_var', 'batch'], capture_output = True, text = True, timeout = 30)
  assert invalid.returncode != 0
  assert 'test_group_var and nonempty test_groups' in invalid.stderr

  # Numerical settings remain mandatory instead of duplicating schema defaults.
  seed_flag = command.index('--seed')
  no_seed = command[:seed_flag] + command[seed_flag + 2:]
  invalid = subprocess.run(no_seed, capture_output = True, text = True, timeout = 30)
  assert invalid.returncode != 0 and '--seed' in invalid.stderr

  classes = tmp_path / 'classes.csv'
  classes.write_text('class\na\nb\n')
  subprocess.run([sys.executable, str(ROOT / 'scripts/train_xgboost.py'), 'specs',
    '--cls_f', str(classes), '--counts_f', str(settings['counts_f']),
    '--predictions_f', str(tmp_path / 'unused.csv.gz'), '--output_dir', str(tmp_path),
    '--xgboost_tag', 'demo'], check = True, capture_output = True, text = True, timeout = 30)
  assert pl.read_csv(tmp_path / 'demo_coarse_evaluation.csv').is_empty()
  assert set(pl.read_csv(tmp_path / 'demo_label_specs.csv')['label']) == {'a', 'b'}


def test_specification_csv_preserves_strings_and_strips_values(tmp_path):
  fine = tmp_path / 'labels.csv'
  fine.write_text('label,colour,coarse_label\n001, #112233 , 01 \n"a,b",#223344,02\n')
  coarse = tmp_path / 'coarse.csv'
  coarse.write_text('coarse_label,coarse_colour\n 02 ,#123456\n01, #654321 \n')

  fine_specs, coarse_specs = read_label_specs(fine, coarse, ['001', 'a,b'])
  assert [row['label'] for row in fine_specs] == ['001', 'a,b']
  assert [row['coarse_label'] for row in fine_specs] == ['01', '02']
  assert [row['coarse_label'] for row in coarse_specs] == ['02', '01']
  assert fine_specs[0]['colour'] == '#112233'
  assert coarse_specs[1]['coarse_colour'] == '#654321'


@pytest.mark.parametrize('coarse_table', [False, True])
def test_specification_csv_rejects_extra_fields(tmp_path, coarse_table):
    fine = tmp_path / 'fine.csv'
    coarse = tmp_path / 'coarse.csv'
    fine.write_text('label,coarse_label\na,A\n')
    if coarse_table:
        coarse.write_text('coarse_label\nA,unexpected\n')
    else:
        fine.write_text('label,coarse_label\na,A,unexpected\n')
    import polars as pl
    with pytest.raises(pl.exceptions.ComputeError):
        read_label_specs(fine, coarse if coarse_table else None)


def test_evaluation_requires_selection_metadata():
    import polars as pl
    from train_xgboost import _evaluation_counts
    predictions_df = pl.DataFrame({'label': ['a'], 'predicted_label': ['a'], 'split': ['train']})
    with pytest.raises(ValueError, match = 'selected column'):
        _evaluation_counts(predictions_df)
