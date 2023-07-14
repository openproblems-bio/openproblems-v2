import yaml
import scanpy as sc

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/integrated_embedding.h5ad',
    'ouput': 'output.h5ad'
}
## VIASH END

with open(meta['config'], 'r', encoding="utf8") as file:
    config = yaml.safe_load(file)

output_type = config["functionality"]["info"]["subtype"]

print('Read input', flush=True)
adata = sc.read_h5ad(par['input'])

print('Run kNN', flush=True)
sc.pp.neighbors(adata, use_rep='X_emb')

print("Store outputs", flush=True)
adata.uns['output_type'] = output_type
adata.write_h5ad(par['output'], compression='gzip')