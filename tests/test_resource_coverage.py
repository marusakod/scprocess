import json
import re
from pathlib import Path

from scripts.scprocess_utils import (
  RESOURCE_CHEMISTRY_GROUPS, _get_chemistry_group, _load_resource_models)


ROOT = Path(__file__).resolve().parents[1]
CONFIG_SCHEMA_F = ROOT / 'resources/schemas/config.schema.json'
JOIN_SCHEMA_F = ROOT / 'resources/schemas/join.schema.json'

# join workflows take their config keys from the join schema, the rest from the project schema
JOIN_RULE_FILES = {'join.smk', 'join_zoom.smk'}

# Shiny data rules were split from the old build_*_shiny_app rules; users' existing
# gb_/mins_ keys for the old names still apply (see _shiny_data_resource in shiny.smk)
SHINY_OLD_RULE_NAMES = {'build_shiny_data': 'build_shiny_app', 'build_zoom_shiny_data': 'build_zoom_shiny_app'}


def _call_arguments(text, start):
  # arguments of the call whose opening parenthesis is at text[start], split at top-level commas
  depth, args, current = 0, [], ''
  for char in text[start:]:
    if char in '([{':
      depth += 1
      if depth == 1:
        continue
    elif char in ')]}':
      depth -= 1
      if depth == 0:
        args.append(current.strip())
        return args
    if char == ',' and depth == 1:
      args.append(current.strip())
      current = ''
    else:
      current += char
  raise ValueError('unbalanced parentheses')


def _find_resource_calls():
  # {(rule, param): {rule file: [argument lists]}} for every get_resources() call and
  # _shiny_data_resource() call (which forwards to get_resources)
  calls = {}
  for rule_f in sorted((ROOT / 'rules').glob('*.smk')):
    text = rule_f.read_text()

    for match in re.finditer(r'get_resources\(', text):
      args = _call_arguments(text, match.end() - 1)
      if args[0] != 'RESOURCE_PARAMS':
        continue
      # args: RESOURCE_PARAMS, rules, input, rule, param, attempt[, run]
      rule, param = args[3].strip("'"), args[4].strip("'")
      calls.setdefault((rule, param), {}).setdefault(rule_f.name, []).append(args)

    for match in re.finditer(r'_shiny_data_resource\(', text):
      args = _call_arguments(text, match.end() - 1)
      if not args[0].startswith("'"):
        continue
      # args: rule, old_rule, input, param, attempt
      rule, param = args[0].strip("'"), args[3].strip("'")
      calls.setdefault((rule, param), {}).setdefault(rule_f.name, []).append(args)

  return calls


def _direct_config_reads():
  reads = set()
  for rule_f in (ROOT / 'rules').glob('*.smk'):
    for prefix, rule in re.findall(r"config\['resources'\]\['(gb|mins)_(\w+)'\]", rule_f.read_text()):
      reads.add((rule, 'memory' if prefix == 'gb' else 'time'))
  return reads


def _schema_resource_properties(schema_f):
  schema = json.loads(schema_f.read_text())
  properties = schema['properties']['resources']['properties']
  return {key: spec for key, spec in properties.items() if key.startswith(('gb_', 'mins_'))}


def _key(rule, param):
  return ('gb_' if param == 'memory' else 'mins_') + rule


def test_call_finder_sees_the_rules():
  calls = _find_resource_calls()
  assert len(calls) > 150
  assert ('run_cellbender', 'memory') in calls
  assert ('build_zoom_shiny_data', 'time') in calls
  assert ('zoom_train_xgboost_prepare', 'memory') in calls


def test_every_called_rule_has_a_model_or_a_floor():
  models, floors = _load_resource_models(ROOT)
  uncovered = [
    (rule, param) for rule, param in _find_resource_calls()
    if param not in models.get(rule, {}) and param not in floors.get(rule, {})]
  assert uncovered == []


def test_every_model_and_floor_belongs_to_a_called_rule():
  models, floors = _load_resource_models(ROOT)
  calls = _find_resource_calls()
  stale_models = [(rule, param) for rule, params in models.items() for param in params
                  if (rule, param) not in calls]
  stale_floors = [(rule, param) for rule, params in floors.items() for param in params
                  if (rule, param) not in calls]
  assert stale_models == []
  assert stale_floors == []


def test_no_rule_reads_resource_values_from_the_config_directly():
  # such rules ignore the models, floors and retries; use get_resources() instead
  assert _direct_config_reads() == set()


def test_resource_keys_in_the_schemas_have_no_defaults():
  # defaults live in resource_floors.yaml, so that a key in the config is always an override
  for schema_f in (CONFIG_SCHEMA_F, JOIN_SCHEMA_F):
    with_default = [key for key, spec in _schema_resource_properties(schema_f).items()
                    if 'default' in spec]
    assert with_default == [], schema_f.name


def test_every_called_rule_can_be_overridden_in_the_config():
  schema_keys = {
    'config': set(_schema_resource_properties(CONFIG_SCHEMA_F)),
    'join': set(_schema_resource_properties(JOIN_SCHEMA_F))}

  missing = []
  for (rule, param), by_file in _find_resource_calls().items():
    for rule_f in by_file:
      kinds = ['config', 'join'] if rule_f == 'shiny.smk' else ['join' if rule_f in JOIN_RULE_FILES else 'config']
      names = {_key(rule, param), _key(SHINY_OLD_RULE_NAMES.get(rule, rule), param)}
      if not any(names & schema_keys[kind] for kind in kinds):
        missing.append((rule, param, rule_f))
  assert missing == []


def test_every_schema_resource_key_belongs_to_a_rule():
  used = {rule for rule, _ in _find_resource_calls()}
  used |= set(SHINY_OLD_RULE_NAMES.values())
  for schema_f in (CONFIG_SCHEMA_F, JOIN_SCHEMA_F):
    orphans = {key for key in _schema_resource_properties(schema_f)
               if key.split('_', 1)[1] not in used}
    assert orphans == set(), schema_f.name


def test_chemistry_dependent_models_are_called_with_a_run():
  models, _ = _load_resource_models(ROOT)
  chemistry_dependent = {
    (rule, param) for rule, params in models.items() for param, chem_models in params.items()
    if 'all' not in chem_models}
  assert chemistry_dependent

  for pair in chemistry_dependent:
    for rule_f, arg_lists in _find_resource_calls()[pair].items():
      for args in arg_lists:
        # RESOURCE_PARAMS, rules, input, rule, param, attempt, run
        assert len(args) == 7, f'{pair} in {rule_f} is called without a run'


def _find_enum(node, name):
  # first enum listed for the property 'name' anywhere in a schema
  if isinstance(node, dict):
    if name in node.get('properties', {}) and 'enum' in node['properties'][name]:
      return node['properties'][name]['enum']
    for value in node.values():
      found = _find_enum(value, name)
      if found:
        return found
  elif isinstance(node, list):
    for value in node:
      found = _find_enum(value, name)
      if found:
        return found
  return None


def test_every_config_chemistry_maps_to_a_model_group():
  chemistries = _find_enum(json.loads(CONFIG_SCHEMA_F.read_text()), 'tenx_chemistry')
  assert 'flexv1' in chemistries and '5v3' in chemistries
  for chemistry in chemistries:
    assert _get_chemistry_group(chemistry) in RESOURCE_CHEMISTRY_GROUPS
