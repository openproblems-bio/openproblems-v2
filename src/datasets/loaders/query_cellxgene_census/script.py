import os
import sys
import cellxgene_census

## VIASH START
par = {
    "input_uri": None,
    "census_version": "stable",
    "species": "homo_sapiens",
    "obs_value_filter": "is_primary_data == True and cell_type_ontology_term_id in ['CL:0000136', 'CL:1000311', 'CL:0002616'] and suspension_type == 'cell'",
    "cell_filter_grouping": ["dataset_id", "tissue", "assay", "disease", "cell_type"],
    "cell_filter_minimum_count": 100,
    "obs_batch": [ "batch" ],
    "obs_batch_separator": "+",
    "output": "output.h5ad",
    "output_compression": "gzip",
}
meta = {"resources_dir": "src/common/helper_functions"}
## VIASH END

sys.path.append(meta["resources_dir"])

from setup_logger import setup_logger
logger = setup_logger()

def connect_census(uri, census_version):
    """
    Connect to CellxGene Census or user-provided TileDBSoma object
    """
    ver = census_version or "stable"
    logger.info("Connecting to CellxGene Census at %s", f"'{uri}'" if uri else f"version '{ver}'")
    return cellxgene_census.open_soma(uri=uri, census_version=ver)

def get_anndata(census_connection, obs_value_filter, species):
    logger.info("Getting gene expression data based on %s query.", obs_value_filter)
    return cellxgene_census.get_anndata(
        census=census_connection, obs_value_filter=obs_value_filter, organism=species
    )

def add_batch_to_obs(adata, par):
    logger.info("Add batch to the AnnData object.")
    if par["obs_batch"]:
        # fetch batch columns from obs
        cols = [adata.obs[key] for key in par["obs_batch"]]
        
        # join cols
        obs_batch = [par["obs_batch_separator"].join(row) for row in zip(*cols)]

        # store in adata
        adata.obs["batch"] = obs_batch

def add_metadata_to_uns(adata, par):
    logger.info("Add metadata to the AnnData object.")
    for key in ["dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism"]:
        adata.uns[key] = par[key]

def cellcensus_cell_filter(adata, cell_filter_grouping, cell_filter_minimum_count):
    t0 = adata.shape
    adata = adata[
        adata.obs.groupby(cell_filter_grouping)["soma_joinid"].transform("count")
        >= cell_filter_minimum_count
    ]
    t1 = adata.shape
    logger.info(
        "Removed %s cells based on %s cell_filter_minimum_count of %s cell_filter_grouping."
        % ((t0[0] - t1[0]), cell_filter_minimum_count, cell_filter_grouping)
    )
    return adata

def write_anndata(adata, path, compression):
    logger.info("Writing AnnData object to '%s'", path)

    adata.write_h5ad(path, compression=compression)

def print_unique(adata, column):
    formatted = "', '".join(adata.obs[column].unique())
    logger.info(f"Unique {column}: ['{formatted}']")

def print_summary(adata):
    logger.info(f"Resulting dataset: {adata}")

    logger.info("Summary of dataset:")
    print_unique(adata, "assay")
    print_unique(adata, "assay_ontology_term_id")
    print_unique(adata, "cell_type")
    print_unique(adata, "cell_type_ontology_term_id")
    print_unique(adata, "dataset_id")
    print_unique(adata, "development_stage")
    print_unique(adata, "development_stage_ontology_term_id")
    print_unique(adata, "disease")
    print_unique(adata, "disease_ontology_term_id")
    print_unique(adata, "tissue")
    print_unique(adata, "tissue_ontology_term_id")
    print_unique(adata, "tissue_general")
    print_unique(adata, "tissue_general_ontology_term_id")

def main():
    # check arguments
    if (par["cell_filter_grouping"] is None) != (par["cell_filter_minimum_count"] is None):
        raise NotImplementedError(
            "You need to specify either both or none of the following parameters: cell_filter_grouping, cell_filter_minimum_count"
        )
    
    with connect_census(uri=par["input_uri"], census_version=par["census_version"]) as conn:
        adata = get_anndata(conn, par["obs_value_filter"], par["species"])

    if par["cell_filter_grouping"] is not None:
        adata = cellcensus_cell_filter(
            adata,
            par["cell_filter_grouping"],
            par["cell_filter_minimum_count"]
        )

    # use feature_id as var_names
    adata.var_names = adata.var["feature_id"]

    add_batch_to_obs(adata, par)

    # add metadata to uns
    add_metadata_to_uns(adata, par)

    # print summary
    print_summary(adata)

    # write output to file
    write_anndata(adata, par["output"], par["output_compression"])


if __name__ == "__main__":
    main()
