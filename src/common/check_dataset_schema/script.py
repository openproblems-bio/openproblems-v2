import anndata as ad
import yaml
import shutil
import json
import numpy as np
import pandas as pd
import scipy

## VIASH START
par = {
  'input': 'resources_test/common/bmmc_multiome_starter/dataset_rna.h5ad',
  'schema': None,
  'stop_on_error': False,
  'checks': None,
  'output': None,
  'meta': 'output/meta.json',
}
## VIASH END

def check_structure(slot_info, adata_slot):
  missing = []
  for obj in slot_info:
    if obj.get('required') and obj['name'] not in adata_slot:
      missing.append(obj['name'])
  return missing

print('Load data', flush=True)
adata = ad.read_h5ad(par['input']).copy()

# create data structure
out = {
  "exit_code": 0,
  "error": {},
  "data_schema": "ok"
}

def is_atomic(obj):
  return isinstance(obj, str) or isinstance(obj, int) or isinstance(obj, bool) or isinstance(obj, float)

def to_atomic(obj):
  if isinstance(obj, np.float64):
    return float(obj)
  elif isinstance(obj, np.int64):
    return int(obj)
  elif isinstance(obj, np.bool_):
    return bool(obj)
  elif isinstance(obj, np.str_):
    return str(obj)
  return obj

def is_list_of_atomics(obj):
  if not isinstance(obj, (list,pd.core.series.Series,np.ndarray)):
    return False
  return all(is_atomic(elem) for elem in obj)

def to_list_of_atomics(obj):
  if isinstance(obj, pd.core.series.Series):
    obj = obj.to_numpy()
  if isinstance(obj, np.ndarray):
    obj = obj.tolist()
  return [to_atomic(elem) for elem in obj]

def is_dict_of_atomics(obj):
  if not isinstance(obj, dict):
    return False
  return all(is_atomic(elem) for _, elem in obj.items())

def to_dict_of_atomics(obj):
  return {k: to_atomic(v) for k, v in obj.items()}

def get_structure_shape(obj) -> list:
  if isinstance(obj, np.ndarray):
    return list(obj.shape)
  elif scipy.sparse.issparse(obj):
    return list(obj.shape)
  elif isinstance(obj, pd.core.frame.DataFrame):
    return list(obj.shape)
  elif isinstance(obj, pd.core.series.Series):
    return list(obj.shape)
  elif isinstance(obj, list):
    return list(len(obj))
  elif isinstance(obj, dict):
    return list(len(obj))
  elif is_atomic(obj):
    return [1]
  return None

def get_structure_type(obj) -> str:
  # return one of: atomic, dataFrame, vector, dict, denseMatrix, sparseMatrix
  if is_atomic(obj):
    return "atomic"
  elif isinstance(obj, (list,pd.core.series.Series)):
    return "vector"
  elif isinstance(obj, dict):
    return "dict"
  elif isinstance(obj, pd.core.frame.DataFrame):
    return "dataframe"
  elif scipy.sparse.issparse(obj):
    return "sparsematrix"
  elif isinstance(obj, np.ndarray):
    return "densematrix"
  return "other: " + str(type(obj))

def get_structure_dtype(obj) -> str:
  if isinstance(obj, np.ndarray):
    return obj.dtype.name
  elif isinstance(obj, pd.core.series.Series):
    return obj.dtype.name
  elif isinstance(obj, pd.core.frame.DataFrame):
    return [dtype.name for dtype in obj.dtypes]
  elif scipy.sparse.issparse(obj):
    return obj.dtype.name
  elif is_atomic(obj):
    return type(obj).__name__
  return None

def get_structure(adata, struct):
  adata_struct = getattr(adata, struct)
  if (struct == "X"):
    adata_struct = {"X": adata_struct} if adata_struct is not None else {}
  return [
    {
      "name": key,
      "type": get_structure_type(value),
      "shape": get_structure_shape(value),
      "dtype": get_structure_dtype(value),
    }
    for key, value in adata_struct.items()
  ]

if par['meta'] is not None:
  print("Extract metadata from object", flush=True)
  uns = {}
  for key, val in adata.uns.items():
    if is_atomic(val):
      uns[key] = to_atomic(val)
    elif is_list_of_atomics(val) and len(val) <= 10:
      uns[key] = to_list_of_atomics(val)
    elif is_dict_of_atomics(val) and len(val) <= 10:
      uns[key] = to_dict_of_atomics(val)
  structure = {
    struct: get_structure(adata, struct)
    for struct
    in ["X", "obs", "var", "obsp", "varp", "obsm", "varm", "layers", "uns"]
  }
  meta = {"uns": uns, "structure": structure}
  with open(par["meta"], "w") as f:
    yaml.dump(meta, f, indent=2)

if par['schema'] is not None:
  print("Check AnnData against schema", flush=True)
  with open(par["schema"], "r") as f:
    data_struct = yaml.safe_load(f)

  def_slots = data_struct['info']['slots']

  missing= []
  for slot in def_slots:
    missing_x = False
    if slot == "X":
      if adata.X is None:
        missing_x = True
      continue
    missing = check_structure(def_slots[slot], getattr(adata, slot))
    if missing_x:
      missing.append("X")
    if missing:
      out['exit_code'] = 1
      out['data_schema'] = 'not ok'
      out['error'][slot] = missing

  if par['checks'] is not None:
    with open(par["checks"], "w") as f:
      json.dump(out, f, indent=2)

if par['output'] is not None and out["data_schema"] == "ok":
  shutil.copyfile(par["input"], par["output"])

if par['stop_on_error']:
  exit(out['exit_code'])  
