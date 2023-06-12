from pprint import pprint
import pathlib
import logging
from jsonschema import validate
import yaml

## VIASH START

meta = {
    "config" : "src/tasks/batch_integration/methods/bbknn/config.vsh.yaml",
    "resources_dir": "src/common/check_yaml_schema"
}

## VIASH END

pathlib.Path(par['output']).parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(par['output'], 'w'),
        logging.StreamHandler()
    ]
)

input_yaml_file = meta["config"]
schema_yaml_file = meta["resources_dir"] + "/schema.yaml"

print('Read files...', flush=True)
with open(input_yaml_file, 'r') as f:
  input_yaml = yaml.safe_load(f)
  print('Input YAML:')
  pprint(input_yaml)

with open(schema_yaml_file, 'r') as f:
  schema_yaml = yaml.safe_load(f)
  print('Schema YAML:')
  pprint(schema_yaml)

print('Validate schema...', flush=True)
try:
  validate(input_yaml, schema_yaml)
except Exception as e:
  logging.exception('YAML validation failed')
  if par['stop_on_error']:
    exit(1)
else:
  print('Validator passed')
