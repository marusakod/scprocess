"""Scientific split, sampling and count-reading regressions for native training."""
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest
import yaml

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from train_xgboost_plan import balanced_sample, assign_splits, prepare_plan
from train_xgboost_workflow import read_selected_h5ad, best_iteration_model, aggregate_plan_predictions


def frame():
 rows=[]
 for donor in ['d1','d2','d3']:
  for tube in range(5):
   for label in ['a','b']:
    for cell in range(10):
     rows.append(dict(cell_id=f'{donor}_{tube}_{label}_{cell}',label=label,donor=donor,
      tube=f'{donor}_t{tube}',project='p',batch=f'{donor}_{tube}'))
 return pd.DataFrame(rows)


def cfg():
 return dict(test_group_var='donor',test_groups=['d3'],validation_group_var='tube',
  validation_fraction=.4,validation_variable_combinations=['donor','project'],min_cells_per_type=5,
  excluded_labels=[],train_cells_per_type=15,validation_cells_per_type=8,
  sampling_variable_combinations=['donor','tube'],seed=7)


def test_donor_and_tube_disjoint():
 result,_,_=assign_splits(frame(),cfg(),np.random.default_rng(7))
 assert set(result[result.split=='test'].donor)=={'d3'}
 assert 'd3' not in set(result[result.split!='test'].donor)
 assert not set(result[result.split=='train'].tube)&set(result[result.split=='validation'].tube)
 assert set(result.split)=={'train','validation','test'}


def test_group_spanning_combinations_fails():
 data=frame();data.loc[data.donor=='d2','tube']='shared'
 data.loc[data.donor=='d1','tube']='shared'
 with pytest.raises(ValueError,match='spans'):assign_splits(data,cfg(),np.random.default_rng(7))


def test_unknown_test_group_fails():
 config=cfg();config['test_groups']=['absent']
 with pytest.raises(ValueError,match='absent'):assign_splits(frame(),config,np.random.default_rng(7))


def test_joint_combination_sampling():
 data=pd.DataFrame({'label':['a']*105,'group':['large']*100+['rare1']*3+['rare2']*2})
 sample=balanced_sample(data,9,['group'],np.random.default_rng(8))
 assert len(sample)==9 and sample.index.is_unique
 assert sample.groupby('group').size().to_dict()=={'large':4,'rare1':3,'rare2':2}


def test_prepare_plan_and_no_default_smoothing(tmp_path):
 data=frame().rename(columns={'label':'annotation'})
 annotations=tmp_path/'annotations.csv';data.to_csv(annotations,index=False)
 integration=tmp_path/'integration.csv.gz'
 data[['cell_id','batch']].rename(columns={'batch':'sample_id'}).assign(**{'RNA_snn_res.2':'cl1'}).to_csv(integration,index=False)
 manifest=tmp_path/'manifest.yaml';manifest.write_text(yaml.safe_dump({b:'unused.h5ad' for b in data.batch.unique()}))
 config=cfg();config['excluded_labels']=['b']
 plan=prepare_plan(annotations,integration,manifest,config,'sample_id',tmp_path/'plan.csv.gz',tmp_path/'coverage.csv',tmp_path/'counts.csv')
 assert set(plan[plan.selected].label)=={'a'}
 assert set(plan[plan.label=='b'].split)=={'excluded'}
 assert plan.loc[plan.label=='a','label'].equals(plan.loc[plan.label=='a','original_label'])
 assert len(plan[plan.selected & (plan.split=='train')])==15
 assert len(plan[plan.selected & (plan.split=='validation')])==8


def test_schema_renames_and_test_group_pair():
 from jsonschema import Draft202012Validator
 for context in ['config','join','zoom']:
  schema=json.loads((ROOT/f'resources/schemas/{context}.schema.json').read_text())['properties']['train_xgboost']
  good=dict(annots_f='annotations.csv',xgboost_tag='model',**cfg())
  # cfg uses internal donor/project metadata names but the same public keys.
  good.pop('excluded_labels');good['excluded_labels']=['Misc neuron']
  Draft202012Validator(schema).validate(good)
  for bad in [dict(annots_f='a',ref_tag='old'),dict(annots_f='a',refine_labels=True),
              dict(annots_f='a',test_group_var='donor')]:
   assert list(Draft202012Validator(schema).iter_errors(bad))


def test_backed_rows_align_and_preserve_bytes(tmp_path):
 import anndata as ad
 import scipy.sparse as sp
 paths=[]
 for name,genes,raw in [('a',['g1','g2'],[[1,2],[3,4]]),('b',['g2','g1'],[[6,5],[8,7]])]:
  path=tmp_path/f'{name}.h5ad'
  data=ad.AnnData(sp.csr_matrix(raw));data.obs_names=['c1','c2'];data.var_names=genes
  data.write_h5ad(path);paths.append(path)
 before=[(p.read_bytes(),p.stat().st_mtime_ns) for p in paths]
 raw,genes,_=read_selected_h5ad(paths[0],['c2','c1'])
 assert raw.toarray().tolist()==[[3,4],[1,2]]
 with pytest.warns(UserWarning,match='different order'):
  raw,_,_=read_selected_h5ad(paths[1],['c1'],genes)
 assert raw.toarray().tolist()==[[5,6]]
 assert before==[(p.read_bytes(),p.stat().st_mtime_ns) for p in paths]
 with pytest.raises(ValueError,match='Missing'):read_selected_h5ad(paths[0],['absent'])
 with pytest.raises(ValueError,match='Different'):read_selected_h5ad(paths[1],['c1'],['g1','absent'])


def test_best_iteration_export():
 import xgboost as xgb
 x=np.arange(80,dtype=float).reshape(40,2);y=(x[:,0]>30).astype(int)
 model=xgb.train({'objective':'binary:logistic','nthread':1,'seed':7},xgb.DMatrix(x,label=y),num_boost_round=6)
 model.set_attr(best_iteration='2')
 trimmed=best_iteration_model(model)
 assert trimmed.num_boosted_rounds()==3
 np.testing.assert_array_equal(trimmed.predict(xgb.DMatrix(x)),model.predict(xgb.DMatrix(x),iteration_range=(0,3)))


def test_aggregation_preserves_test_and_counts_unsupported(tmp_path):
 plan=tmp_path/'plan.csv.gz'
 pd.DataFrame(dict(cell_id=['a','b','c'],label=['x','unseen','x'],original_label=['x','unseen','x'],
  split=['train','test','validation'],selected=[True,False,False],label_supported_in_training=[True,False,True])).to_csv(plan,index=False)
 pred=tmp_path/'pred.csv.gz'
 pd.DataFrame(dict(cell_id=['a','b','c','outside'],predicted_label_naive=['x','x','x','x'],probability_naive=[.9,.8,.9,.7])).to_csv(pred,index=False)
 integration=tmp_path/'integration.csv';pd.DataFrame({'cell_id':['a','b','c']}).to_csv(integration,index=False)
 output=tmp_path/'demo_fulldata_predictions.csv.gz'
 aggregate_plan_predictions([pred],plan,output,integration)
 result=pd.read_csv(output)
 assert set(result.cell_id)=={'a','b','c'}
 assert result.set_index('cell_id').loc['b','split']=='test'
 metrics=pd.read_csv(tmp_path/'demo_evaluation.csv').set_index('split')
 assert metrics.loc['test','accuracy']==0 and metrics.loc['test','macro_f1']==0
 assert 'development_holdout' in metrics.index


def test_validation_only_genes_do_not_pass_detection_filter():
 from train_xgboost_workflow import training_gene_mask
 import scipy.sparse as sp
 raw=sp.csr_matrix([[1,0,1],[2,0,0],[3,0,0],[0,4,0],[0,5,0]])
 mask=training_gene_mask(raw,np.array([True,True,True,False,False]),2)
 assert mask.tolist()==[True,False,False]


def test_majority_votes_use_training_labels_only(tmp_path):
 data=frame().rename(columns={'label':'annotation'})
 data.loc[data.donor!='d3','annotation']='a'
 minority=(data.donor!='d3')&data.cell_id.str.endswith('_b_9')
 data.loc[minority,'annotation']='b';data.loc[data.donor=='d3','annotation']='b'
 annotations=tmp_path/'annotations.csv';data.to_csv(annotations,index=False)
 integration=tmp_path/'integration.csv'
 data[['cell_id','batch']].rename(columns={'batch':'sample_id'}).assign(**{'RNA_snn_res.2':'cl1'}).to_csv(integration,index=False)
 manifest=tmp_path/'manifest.yaml';manifest.write_text(yaml.safe_dump({b:'unused.h5ad' for b in data.batch.unique()}))
 config=cfg();config.update(cluster_majority_vote_resolution=2,cluster_majority_vote_min_purity=.9)
 plan=prepare_plan(annotations,integration,manifest,config,'sample_id',tmp_path/'plan.csv.gz',tmp_path/'coverage.csv',tmp_path/'counts.csv')
 assert set(plan[plan.split=='train'].label)=={'a'}
 assert set(plan[plan.split=='test'].label)=={'b'}
 assert set(plan[plan.split=='test'].original_label)=={'b'}


def test_plan_configuration_keeps_label_map_but_not_hyperparameters():
 from scprocess_utils import get_xgboost_planning_config
 config={'label_map_f':'coarse.csv','seed':3,'pass1_nrounds':12}
 assert get_xgboost_planning_config(config)=={'label_map_f':'coarse.csv','seed':3}


def test_gpu_request_rejects_cpu_build(monkeypatch):
 import train_xgboost_workflow as workflow
 monkeypatch.setattr(workflow.xgb,'build_info',lambda:{'USE_CUDA':False})
 with pytest.raises(RuntimeError,match='without CUDA'):workflow.verify_gpu(2)
