import numpy as np
import anndata as ad
import json
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.model import TransformerModel
from scgpt.utils.util import load_pretrained
import torch

## VIASH START

par = {
    'input': 'resources_test/batch_integration/pancreas/dataset.h5ad',
    'model': 'resources_test/batch_integration/scgpt/pretrained_model',
    'model_config': 'resources_test/batch_integration/scgpt/pretrained_model/config.json',
    'model_vocab': 'resources_test/batch_integration/scgpt/pretrained_model/vocab.json',
    'output': 'output.h5ad'
}

meta = {
    'functionality_name' : 'scgpt',
}

## VIASH END


