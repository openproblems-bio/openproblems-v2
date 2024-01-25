import sys
import anndata as ad
import numpy as np
from pathlib import Path
import torch
import pytorch_lightning as pl
from torch.utils.data import TensorDataset,DataLoader
from pytorch_lightning.callbacks import ModelCheckpoint
from pytorch_lightning.loggers import TensorBoardLogger,WandbLogger

## VIASH START
par = {
    'input_train_mod1': 'output/datasets/predict_modality/openproblems_bmmc_cite_phase2_rna/openproblems_bmmc_cite_phase2_rna.censor_dataset.output_train_mod1.h5ad',
    'input_train_mod2': 'output/datasets/predict_modality/openproblems_bmmc_cite_phase2_rna/openproblems_bmmc_cite_phase2_rna.censor_dataset.output_train_mod2.h5ad',
    'output': 'output/model.pt'
}

meta = {
    'resources_dir': '/helper_functions',
}
## VIASH END


sys.path.append(meta['resources_dir'])
from models import MLP
import utils

def _train(X, y, Xt, yt, enable_ckpt, logger, yaml_path):
  config = utils.load_yaml(yaml_path)
  X = torch.from_numpy(X).float()
  y = torch.from_numpy(y).float()
  ymean = torch.mean(y,dim=0,keepdim=True)
  
  tr_ds = TensorDataset(X,y)
  nw = 0 if 'ATAC' in yaml_path else 4 
  tr_loader = DataLoader(tr_ds, batch_size=config.batch_size,num_workers=nw,
                      shuffle=True, drop_last=True)
  
  Xt = torch.from_numpy(Xt).float()
  yt = torch.from_numpy(yt).float()
  te_ds = TensorDataset(Xt,yt)
  te_loader = DataLoader(te_ds, batch_size=config.batch_size,num_workers=0,
                      shuffle=False, drop_last=False)
  
  checkpoint_callback = ModelCheckpoint(monitor='valid_RMSE')
  if enable_ckpt:
      epochs = config.epochs
      cb = [checkpoint_callback]
  else:
      epochs = 1
      cb = None
  
  trainer = pl.Trainer(enable_checkpointing=enable_ckpt, logger=logger, 
                        gpus=1, max_epochs=epochs, 
                        callbacks=cb,
                        progress_bar_refresh_rate=5)
  
  net = MLP(X.shape[1],y.shape[1],ymean,config)
  trainer.fit(net, tr_loader, te_loader)
  
  cp = 'best' if enable_ckpt else None
  yp = trainer.predict(net,te_loader,ckpt_path=cp)
  yp = torch.cat(yp,dim=0)
  
  score = ((yp-yt)**2).mean()**0.5
  print(f"VALID RMSE {score:.3f}")
  del trainer
  return score,yp.detach().numpy()


def train(task,cp,wp,tr1,tr2):
  yaml_path = f'{cp}/yaml/mlp_{task}.yaml'
  yps = []
  scores = []

  msgs = {}
  for fold in range(3):

    run_name = f"{task}_fold_{fold}"
    save_path = f'{wp}/{run_name}'
    Path(save_path).mkdir(parents=True, exist_ok=True)   

    X,y,Xt,yt = utils.split(tr1, tr2, fold)
    run_name = f'fold_{fold}'
    logger = TensorBoardLogger(save_path, name='') 
    
    enable_ckpt = True
    
    score, yp = _train(X, y, Xt, yt, enable_ckpt, logger, yaml_path)
    yps.append(yp)
    scores.append(score)
    msg = f"{task} Fold {fold} RMSE {score:.3f}"
    msgs[f'Fold {fold}'] = f'{score:.3f}'
    print(msg)

  yp = np.concatenate(yps)
  score = np.mean(scores)
  msgs['Overall'] = f'{score:.3f}'
  print('Overall', f'{score:.3f}')


print('Loading data', flush=True)
input_train_mod1 = ad.read_h5ad(par['input_train_mod1'])
input_train_mod2 = ad.read_h5ad(par['input_train_mod2'])

mod_1 = input_train_mod1.uns['modality']
mod_2 = input_train_mod2.uns['modality']

task = f'{mod_1}2{mod_2}'
train(task,cp=meta['resources_dir'],
      wp=par['output_pretrain'],
      tr1=input_train_mod1,
      tr2=input_train_mod2)