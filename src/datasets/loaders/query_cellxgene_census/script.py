import cellxgene_census

## VIASH START
par = {
    "input_database": "CellxGene",
    "cellxgene_release": "2023-07-25",
    "species": "homo_sapiens",
    "cell_query": "is_primary_data == True and cell_type_ontology_term_id in ['CL:0000136', 'CL:1000311', 'CL:0002616'] and suspension_type == 'cell'",
    "output": "output.h5ad",
}
meta = {}
### VIASH END


def connect_census(input_database, release):
    """
    Connect to CellxGene Census or user-provided TileDBSoma object
    """
    if input_database != "CellxGene":
        raise NotImplementedError("Custom census database is not implemented yet!")
    return cellxgene_census.open_soma(census_version=release)


def get_anndata(census_connection, cell_query, species):
    return cellxgene_census.get_anndata(
        census=census_connection, obs_value_filter=cell_query, organism=species
    )


def add_cellcensus_metadata_obs(census_connection, query_data):
    census_datasets = (
        census_connection["census_info"]["datasets"].read().concat().to_pandas()
    )

    query_data.obs.dataset_id = query_data.obs.dataset_id.astype("category")

    dataset_info = (
        census_datasets[
            census_datasets.dataset_id.isin(query_data.obs.dataset_id.cat.categories)
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

    return query_data.obs.merge(dataset_info, on="dataset_id", how="left")


def cellcensus_cell_filter(query_data, cells_filter_columns, min_cells_filter_columns):
    t0 = query_data.shape
    query_data = query_data[
        query_data.obs.groupby(cells_filter_columns)["soma_joinid"].transform("count")
        >= min_cells_filter_columns
    ]
    t1 = query_data.shape
    return query_data


census_connection = connect_census(par["input_database"], par["cellxgene_release"])

query_data = get_anndata(census_connection, par["cell_query"], par["species"])

query_data.obs = add_cellcensus_metadata_obs(census_connection, query_data)

census_connection.close()
del census_connection

if par["cells_filter_columns"]:
    if not par["min_cells_filter_columns"]:
        raise NotImplementedError(
            "You specified cells_filter_columns, thus add min_cells_filter_columns!"
        )
    query_data = cellcensus_cell_filter(
        query_data, par["cells_filter_columns"], par["min_cells_filter_columns"]
    )

query_data.var_names = query_data.var["feature_id"]
query_data.var["gene_symbol"] = query_data.var["feature_name"]

with open(par["output"], "w") as f:
    query_data.write_h5ad(par["output"], compression="gzip")
