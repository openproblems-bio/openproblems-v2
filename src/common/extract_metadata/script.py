import anndata as ad
import yaml
import numpy as np
import pandas as pd
import scipy
import os
import datetime

## VIASH START
par = {
  'input': 'resources_test/common/pancreas/dataset.h5ad',
  'output': 'output/meta.yaml'
}
## VIASH END

# Get file size in bytes
file_size = os.path.getsize(par['input'])

# Get file creation time
creation_time = os.path.getctime(par['input'])
# Convert creation time from seconds since epoch to a readable timestamp
creation_time = datetime.datetime.fromtimestamp(creation_time)
# Format the datetime object as 'DD-MM-YYYY'
creation_time = creation_time.strftime('%d-%m-%Y')

print('Load data', flush=True)
adata = ad.read_h5ad(par['input']).copy()

# create data structure
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
    return [len(obj)]
  elif isinstance(obj, dict):
    return [len(obj)]
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


print("Extract metadata from object", flush=True)
uns = {}
for key, val in adata.uns.items():
  if is_atomic(val):
    uns[key] = to_atomic(val)
  elif is_list_of_atomics(val) and len(val) <= 10:
    uns[key] = to_list_of_atomics(val)
  elif is_dict_of_atomics(val) and len(val) <= 10:
    uns[key] = to_dict_of_atomics(val)

uns["file_size"] = file_size
uns["date_created"] = str(creation_time)
structure = {
  struct: get_structure(adata, struct)
  for struct
  in ["X", "obs", "var", "obsp", "varp", "obsm", "varm", "layers", "uns"]
}
meta = {"uns": uns, "structure": structure}

print("Write metadata to file", flush=True)
with open(par["output"], "w") as f:
  yaml.dump(meta, f, indent=2)
