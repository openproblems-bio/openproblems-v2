import cellxgene_census

## VIASH START

par = {
    'dataset_id': 'placeholder',
    "cellxgene_release": "2023-07-25",
    'ouput': 'output.h5ad'
}

meta = {
}

## VIASH END

cellxgene_census.download_source_h5ad(par['dataset_id'], to_path=par['output'], census_version=par['cellxgene_release'])
