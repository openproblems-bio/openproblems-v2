import anndata as ad
import sys
from scipy.sparse import csc_matrix
import numpy as np

## VIASH START
par = {
    'input_train_mod1': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/swap/train_mod1.h5ad',
    'input_train_mod2': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/swap/train_mod2.h5ad',
    'input_test_mod1': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/swap/test_mod1.h5ad',
    'input_model': 'output/model',
    'output': 'output/prediction'
}
meta = {
    'resources_dir': 'src/tasks/predict_modality/methods/simple_mlp',
    'cpus': 10
}
## VIASH END

resources_dir = f"{meta['resources_dir']}/resources"
sys.path.append(resources_dir)
from predict import predict

def get_y_dim(task):
  if task == "ADT2GEX":
      return 13953,
  elif task == "GEX2ADT":
    return 134
  elif task == "ATAC2GEX":
    return 13431
  elif task == "GEX2ATAC":
    return 10000

print('Load data', flush=True)
input_train_mod1 = ad.read_h5ad(par['input_train_mod1'])
input_train_mod2 = ad.read_h5ad(par['input_train_mod2'])
input_test_mod1 = ad.read_h5ad(par['input_test_mod1'])

# determine variables
mod_1 = input_train_mod1.uns['modality']
mod_2 = input_train_mod2.uns['modality']

task = f'{mod_1}2{mod_2}'

y_dim = input_test_mod1.shape[1]

ymean = np.asarray(input_train_mod2.layers["normalized"].mean(axis=0))

print('Start predict', flush=True)
if task == 'GEX2ATAC':
    y_pred = ymean*np.ones([input_test_mod1.shape[0],y_dim])
else:
    y_pred = predict(
       ymean=ymean,
       ad=input_test_mod1.layers["normalized"].toarray(),
       task=task,
       y_dim=y_dim,
       folds=[0,1,2],
       cp=resources_dir,
       wp=par['input_model']
    )

y_pred = csc_matrix(y_pred)

adata = ad.AnnData(
    layers=list(normalized=y_pred),
    shape=y_pred.shape,
    uns={
        'dataset_id': input_train_mod1.uns['dataset_id'],
        'method_id': meta['functionality_name'],
    },
)

print('Write data', flush=True)
adata.write_h5ad(par['output'], compression = "gzip")