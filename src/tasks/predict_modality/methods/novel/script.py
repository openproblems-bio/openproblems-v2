import os
import sys

import torch
from torch.utils.data import DataLoader
from sklearn.model_selection import train_test_split

import anndata as ad
import pickle
import numpy as np
from scipy.sparse import csc_matrix

#check gpu available
if (torch.cuda.is_available()):
    device = 'cuda:0' #switch to current device
    print('current device: gpu', flush=True)
else:
    device = 'cpu'
    print('current device: cpu', flush=True)


## VIASH START

par = {
    'input_train_mod1': 'resources_test/predict_modality/neurips2021_bmmc_cite/train_mod1.h5ad',
    'input_train_mod2': 'resources_test/predict_modality/neurips2021_bmmc_cite/train_mod2.h5ad',
    'input_test_mod1': 'resources_test/predict_modality/neurips2021_bmmc_cite/test_mod1.h5ad',
    'input_test_mod2': 'resources_test/predict_modality/neurips2021_bmmc_cite/test_mod2.h5ad',
    'pretrain': 'model.pt'
}
meta = {
    'resources_dir': '.',
    'functionality_name': '171129'
}
## VIASH END

sys.path.append(meta['resources_dir'])
from helper_functions import train_and_valid, lsiTransformer, ModalityMatchingDataset
from helper_functions import ModelRegressionAtac2Gex, ModelRegressionAdt2Gex, ModelRegressionGex2Adt, ModelRegressionGex2Atac

print("Load data", flush=True)

input_train_mod1 = ad.read_h5ad(par['input_train_mod1'])
input_train_mod2 = ad.read_h5ad(par['input_train_mod2'])
input_test_mod1 = ad.read_h5ad(par['input_test_mod1'])
input_test_mod2 = ad.read_h5ad(par['input_test_mod2'])

mod1 = input_train_mod1.var['feature_types'][0]
mod2 = input_train_mod2.var['feature_types'][0]

model_fp = meta['resources_dir'] + '/' + par['pretrain']

print("Start train", flush=True)
if mod1 != "ADT":
    input_train_mod2_df = input_train_mod2.to_df()
    
    lsi_transformer_gex = lsiTransformer(n_components=256)
    gex_train = lsi_transformer_gex.fit_transform(input_train_mod1)
    
    train_mod1, test_mod1, train_mod2, test_mod2 = train_test_split(gex_train, input_train_mod2_df, test_size=0.25, random_state=666)
    input_train_mod2_df = input_train_mod2.to_df()
else:
    train_mod1 = input_train_mod1.to_df()
    train_mod2 = input_train_mod2.to_df()
    test_mod1 = input_test_mod1.to_df()
    test_mod2 = input_test_mod2.to_df()


if mod1 == 'ATAC' and mod2 == 'GEX':
    dataset_train = ModalityMatchingDataset(train_mod1, train_mod2)
    dataloader_train = DataLoader(dataset_train, 256, shuffle = True, num_workers = 8)

    dataset_test = ModalityMatchingDataset(test_mod1, test_mod2)
    dataloader_test = DataLoader(dataset_test, 64, shuffle = False, num_workers = 8)

    model = ModelRegressionAtac2Gex(256,13431).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=0.00008386597445284492,weight_decay=0.000684887347727808)
        
elif mod1 == 'ADT' and mod2 == 'GEX':
    dataset_train = ModalityMatchingDataset(train_mod1, train_mod2)
    dataloader_train = DataLoader(dataset_train, 64, shuffle = True, num_workers = 4)

    dataset_test = ModalityMatchingDataset(test_mod1, test_mod2)
    dataloader_test = DataLoader(dataset_test, 32, shuffle = False, num_workers = 4)

    model = ModelRegressionAdt2Gex(134,13953).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.00041, weight_decay=0.0000139)


elif mod1 == 'GEX' and mod2 == 'ADT':
    dataset_train = ModalityMatchingDataset(train_mod1, train_mod2)
    dataloader_train = DataLoader(dataset_train, 32, shuffle = True, num_workers = 8)

    dataset_test = ModalityMatchingDataset(test_mod1, test_mod2)
    dataloader_test = DataLoader(dataset_test, 64, shuffle = False, num_workers = 8)

    model = ModelRegressionGex2Adt(256,134).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=0.000034609210829678734, weight_decay=0.0009965881574697426)


elif mod1 == 'GEX' and mod2 == 'ATAC':
    dataset_train = ModalityMatchingDataset(train_mod1, train_mod2)
    dataloader_train = DataLoader(dataset_train, 64, shuffle = True, num_workers = 8)

    dataset_test = ModalityMatchingDataset(test_mod1, test_mod2)
    dataloader_test = DataLoader(dataset_test, 64, shuffle = False, num_workers = 8)

    model = ModelRegressionGex2Atac(256,10000).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=0.00001806762345275399, weight_decay=0.0004084171379280058)

loss_fn = torch.nn.MSELoss()
train_and_valid(model, optimizer, loss_fn, dataloader_train, dataloader_test, model_fp, device)

print("Start predict", flush=True)

if mod1 == 'GEX' and mod2 == 'ADT':
    model = ModelRegressionGex2Adt(256,134)   
    weight = torch.load(model_fp, map_location='cpu')    
        
    model.load_state_dict(weight)    
    input_test_mod1_ = lsi_transformer_gex.transform(input_test_mod1)

elif mod1 == 'GEX' and mod2 == 'ATAC':
    model = ModelRegressionGex2Atac(256,10000)   
    weight = torch.load(model_fp, map_location='cpu')
    with open(par['pretrain'] + '/lsi_transformer.pickle', 'rb') as f:
        lsi_transformer_gex = pickle.load(f)
    
        
    model.load_state_dict(weight)    
    input_test_mod1_ = lsi_transformer_gex.transform(input_test_mod1)
    
elif mod1 == 'ATAC' and mod2 == 'GEX':
    model = ModelRegressionAtac2Gex(256,13431)   
    weight = torch.load(model_fp, map_location='cpu')
    with open(par['pretrain'] + '/lsi_transformer.pickle', 'rb') as f:
        lsi_transformer_gex = pickle.load(f)
        
    model.load_state_dict(weight)    
    input_test_mod1_ = lsi_transformer_gex.transform(input_test_mod1)

elif mod1 == 'ADT' and mod2 == 'GEX':
    model = ModelRegressionAdt2Gex(134,13953)   
    weight = torch.load(model_fp, map_location='cpu')
        
    model.load_state_dict(weight)    
    #input_test_mod1_ = lsi_transformer_gex.transform(input_test_mod1)
    input_test_mod1_ = input_test_mod1.to_df()
    
dataset_test = ModalityMatchingDataset(input_test_mod1_, None, is_train=False)
dataloader_test = DataLoader(dataset_test, 32, shuffle = False, num_workers = 4)

outputs = []
model.eval()
with torch.no_grad():
    for x in dataloader_test:
        output = model(x.float())
        outputs.append(output.detach().cpu().numpy())

outputs = np.concatenate(outputs)
outputs[outputs<0] = 0
outputs = csc_matrix(outputs)

adata = ad.AnnData(
    layers=list(normalized=outputs),
    shape=dim(outputs),
    uns={
        'dataset_id': input_train_mod1.uns['dataset_id'],
        'method_id': meta['functionality_name'],
    },
)
adata.write_h5ad(par['output'], compression = "gzip")


