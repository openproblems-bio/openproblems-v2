import anndata as ad
import sys
from scipy.sparse import csc_matrix
import numpy as np

## VIASH START
par = {
    'input_train_mod1': 'output/datasets/predict_modality/openproblems_bmmc_cite_phase2_rna/openproblems_bmmc_cite_phase2_rna.censor_dataset.output_train_mod1.h5ad',
    'input_train_mod2': 'output/datasets/predict_modality/openproblems_bmmc_cite_phase2_rna/openproblems_bmmc_cite_phase2_rna.censor_dataset.output_train_mod2.h5ad',
    'input_test_mod1': 'output/datasets/predict_modality/openproblems_bmmc_cite_phase2_rna/openproblems_bmmc_cite_phase2_rna.censor_dataset.output_test_mod1.h5ad',
    'input_model': 'path/to/model',
    'output': 'output.h5ad'
}
meta = {
    'resources_dir': 'src/tasks/predict_modality/methods/simplemlp/helper_functions'
}
## VIASH END
sys.path.append(meta['resources_dir'])
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

# Update following becuase it will not work with current setup:
# y_dim,task = get_y_dim(par['input_test_mod1'])
mod_1 = input_train_mod1.uns['modality']
mod_2 = input_train_mod2.uns['modality']

task = f'{mod_1}2{mod_2}'

y_dim = input_test_mod1.shape[1]

ymean = np.asarray(input_train_mod2.layers["counts"].mean(axis=0))

print('Start predict', flush=True)
if task == 'GEX2ATAC':
    y_pred = ymean*np.ones([input_test_mod1.shape[0],y_dim])
else:
    y_pred = predict(ymean,input_test_mod1.layers["counts"].toarray(),task,y_dim,
                     folds=[0,1,2],cp=meta['resources_dir'],
                     wp=par['input_pretrain'])

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