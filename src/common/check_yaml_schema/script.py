import sys
from jsonschema import validate
import yaml

## VIASH START
par = {
  'input': 'src/tasks/batch_integration/methods/bbknn/config.vsh.yaml',
  'schema': 'src/common/api/schema_method.yaml'
}
meta = {
  'functionality_name': 'foo',
}

## VIASH END

input_yaml_file = par.get('input') or meta['config']
schema_yaml_file = par.get('schema') or f"{meta['resources_dir']}/schema.yaml"

with open(input_yaml_file, 'r') as f:
  input_yaml = yaml.safe_load(f)

with open(schema_yaml_file, 'r') as f:
  schema_yaml = yaml.safe_load(f)

validate(input_yaml, schema_yaml)
# try:
#   validate(input_yaml, schema_yaml)
# except Exception as e:
#   sys.exit('YAML validation failed')
