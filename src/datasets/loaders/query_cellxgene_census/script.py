import cellxgene_census

## VIASH START
par = {
    "cellxgene_release": "2023-07-25",
    "species": "homo_sapiens",
    "cell_query": "is_primary_data == True and cell_type_ontology_term_id in ['CL:0000136', 'CL:1000311', 'CL:0002616'] and suspension_type == 'cell'",
    "output": "output.h5ad",
}
meta = {}
### VIASH END

with cellxgene_census.open_soma(census_version=par["cellxgene_release"]) as conn:

    # fetch counts data using query
    query_data = cellxgene_census.get_anndata(
        census=conn,
        obs_value_filter=par["cell_query"],
        organism=par["species"]
    )

    # turn into categorical
    query_data.obs.dataset_id = query_data.obs.dataset_id.astype("category")

    # fetch dataset metadata
    # TODO: fetch only info from the dataset ids in query_data
    dataset_info = (
        conn["census_info"]["datasets"].read().concat().to_pandas()
    )

# subset dataset info
dataset_info = (
    dataset_info[
        dataset_info.dataset_id.isin(query_data.obs.dataset_id.cat.categories)
    ][
        [
            "collection_id",
            "collection_name",
            "collection_doi",
            "dataset_id",
            "dataset_title",
        ]
    ]
    .reset_index(drop=True)
    .apply(lambda x: x.astype("category"))
)

# join dataset info to query data
query_data.obs = query_data.obs.merge(dataset_info, on="dataset_id", how="left")

# use feature id as the var_names
query_data.var_names = query_data.var["feature_id"]

# write to file
with open(par["output"], "w") as f:
    query_data.write_h5ad(par["output"], compression="gzip")
