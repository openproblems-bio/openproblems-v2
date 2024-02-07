# Spatial decomposition

Estimation of cell type proportions per spot in 2D space from spatial
transcriptomic data coupled with corresponding single-cell data

Path:
[`src/tasks/spatial_decomposition`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/tasks/spatial_decomposition)

## Motivation

Spatial decomposition (also often referred to as Spatial deconvolution)
is applicable to spatial transcriptomics data where the transcription
profile of each capture location (spot, voxel, bead, etc.) do not share
a bijective relationship with the cells in the tissue, i.e., multiple
cells may contribute to the same capture location. The task of spatial
decomposition then refers to estimating the composition of cell
types/states that are present at each capture location. The cell
type/states estimates are presented as proportion values, representing
the proportion of the cells at each capture location that belong to a
given cell type.

## Description

We distinguish between *reference-based* decomposition and *de novo*
decomposition, where the former leverage external data (e.g., scRNA-seq
or scNuc-seq) to guide the inference process, while the latter only work
with the spatial data. We require that all datasets have an associated
reference single cell data set, but methods are free to ignore this
information.

## Authors & contributors

| name           | roles              |
|:---------------|:-------------------|
| Giovanni Palla | author, maintainer |
| Scott Gigante  | maintainer         |

## API

``` mermaid
flowchart LR
  file_common_dataset("Common dataset")
  comp_process_dataset[/"Data processor"/]
  file_single_cell("Single cell data")
  file_spatial_masked("Masked")
  file_solution("Solution")
  comp_control_method[/"Control method"/]
  comp_method[/"Method"/]
  comp_metric[/"Metric"/]
  file_output("Output")
  file_score("Score")
  file_common_dataset---comp_process_dataset
  comp_process_dataset-->file_single_cell
  comp_process_dataset-->file_spatial_masked
  comp_process_dataset-->file_solution
  file_single_cell---comp_control_method
  file_single_cell---comp_method
  file_spatial_masked---comp_control_method
  file_spatial_masked---comp_method
  file_solution---comp_control_method
  file_solution---comp_metric
  comp_control_method-->file_output
  comp_method-->file_output
  comp_metric-->file_score
  file_output---comp_metric
```

## File format: Common dataset

A dataset processed by the common dataset processing pipeline.

Example file: `resources_test/common/pancreas/dataset.h5ad`

Description:

This dataset contains both raw counts and normalized data matrices, as
well as a PCA embedding, HVG selection and a kNN graph.

Format:

<div class="small">

    AnnData object
     obs: 'dataset_id', 'assay', 'assay_ontology_term_id', 'cell_type', 'cell_type_ontology_term_id', 'development_stage', 'development_stage_ontology_term_id', 'disease', 'disease_ontology_term_id', 'donor_id', 'is_primary_data', 'organism', 'organism_ontology_term_id', 'self_reported_ethnicity', 'self_reported_ethnicity_ontology_term_id', 'sex', 'sex_ontology_term_id', 'suspension_type', 'tissue', 'tissue_ontology_term_id', 'tissue_general', 'tissue_general_ontology_term_id', 'batch', 'soma_joinid', 'size_factors'
     var: 'feature_id', 'feature_name', 'soma_joinid', 'hvg', 'hvg_score'
     obsm: 'X_pca'
     obsp: 'knn_distances', 'knn_connectivities'
     varm: 'pca_loadings'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism', 'normalization_id', 'pca_variance', 'knn'

</div>

Slot description:

<div class="small">

| Slot                                              | Type      | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
|:--------------------------------------------------|:----------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `obs["dataset_id"]`                               | `string`  | (*Optional*) Identifier for the dataset from which the cell data is derived, useful for tracking and referencing purposes.                                                                                                                                                                                                                                                                                                                                                                    |
| `obs["assay"]`                                    | `string`  | (*Optional*) Type of assay used to generate the cell data, indicating the methodology or technique employed.                                                                                                                                                                                                                                                                                                                                                                                  |
| `obs["assay_ontology_term_id"]`                   | `string`  | (*Optional*) Experimental Factor Ontology (`EFO:`) term identifier for the assay, providing a standardized reference to the assay type.                                                                                                                                                                                                                                                                                                                                                       |
| `obs["cell_type"]`                                | `string`  | (*Optional*) Classification of the cell type based on its characteristics and function within the tissue or organism.                                                                                                                                                                                                                                                                                                                                                                         |
| `obs["cell_type_ontology_term_id"]`               | `string`  | (*Optional*) Cell Ontology (`CL:`) term identifier for the cell type, offering a standardized reference to the specific cell classification.                                                                                                                                                                                                                                                                                                                                                  |
| `obs["development_stage"]`                        | `string`  | (*Optional*) Stage of development of the organism or tissue from which the cell is derived, indicating its maturity or developmental phase.                                                                                                                                                                                                                                                                                                                                                   |
| `obs["development_stage_ontology_term_id"]`       | `string`  | (*Optional*) Ontology term identifier for the developmental stage, providing a standardized reference to the organism’s developmental phase. If the organism is human (`organism_ontology_term_id == 'NCBITaxon:9606'`), then the Human Developmental Stages (`HsapDv:`) ontology is used. If the organism is mouse (`organism_ontology_term_id == 'NCBITaxon:10090'`), then the Mouse Developmental Stages (`MmusDv:`) ontology is used. Otherwise, the Uberon (`UBERON:`) ontology is used. |
| `obs["disease"]`                                  | `string`  | (*Optional*) Information on any disease or pathological condition associated with the cell or donor.                                                                                                                                                                                                                                                                                                                                                                                          |
| `obs["disease_ontology_term_id"]`                 | `string`  | (*Optional*) Ontology term identifier for the disease, enabling standardized disease classification and referencing. Must be a term from the Mondo Disease Ontology (`MONDO:`) ontology term, or `PATO:0000461` from the Phenotype And Trait Ontology (`PATO:`).                                                                                                                                                                                                                              |
| `obs["donor_id"]`                                 | `string`  | (*Optional*) Identifier for the donor from whom the cell sample is obtained.                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `obs["is_primary_data"]`                          | `boolean` | (*Optional*) Indicates whether the data is primary (directly obtained from experiments) or has been computationally derived from other primary data.                                                                                                                                                                                                                                                                                                                                          |
| `obs["organism"]`                                 | `string`  | (*Optional*) Organism from which the cell sample is obtained.                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `obs["organism_ontology_term_id"]`                | `string`  | (*Optional*) Ontology term identifier for the organism, providing a standardized reference for the organism. Must be a term from the NCBI Taxonomy Ontology (`NCBITaxon:`) which is a child of `NCBITaxon:33208`.                                                                                                                                                                                                                                                                             |
| `obs["self_reported_ethnicity"]`                  | `string`  | (*Optional*) Ethnicity of the donor as self-reported, relevant for studies considering genetic diversity and population-specific traits.                                                                                                                                                                                                                                                                                                                                                      |
| `obs["self_reported_ethnicity_ontology_term_id"]` | `string`  | (*Optional*) Ontology term identifier for the self-reported ethnicity, providing a standardized reference for ethnic classifications. If the organism is human (`organism_ontology_term_id == 'NCBITaxon:9606'`), then the Human Ancestry Ontology (`HANCESTRO:`) is used.                                                                                                                                                                                                                    |
| `obs["sex"]`                                      | `string`  | (*Optional*) Biological sex of the donor or source organism, crucial for studies involving sex-specific traits or conditions.                                                                                                                                                                                                                                                                                                                                                                 |
| `obs["sex_ontology_term_id"]`                     | `string`  | (*Optional*) Ontology term identifier for the biological sex, ensuring standardized classification of sex. Only `PATO:0000383`, `PATO:0000384` and `PATO:0001340` are allowed.                                                                                                                                                                                                                                                                                                                |
| `obs["suspension_type"]`                          | `string`  | (*Optional*) Type of suspension or medium in which the cells were stored or processed, important for understanding cell handling and conditions.                                                                                                                                                                                                                                                                                                                                              |
| `obs["tissue"]`                                   | `string`  | (*Optional*) Specific tissue from which the cells were derived, key for context and specificity in cell studies.                                                                                                                                                                                                                                                                                                                                                                              |
| `obs["tissue_ontology_term_id"]`                  | `string`  | (*Optional*) Ontology term identifier for the tissue, providing a standardized reference for the tissue type. For organoid or tissue samples, the Uber-anatomy ontology (`UBERON:`) is used. The term ids must be a child term of `UBERON:0001062` (anatomical entity). For cell cultures, the Cell Ontology (`CL:`) is used. The term ids cannot be `CL:0000255`, `CL:0000257` or `CL:0000548`.                                                                                              |
| `obs["tissue_general"]`                           | `string`  | (*Optional*) General category or classification of the tissue, useful for broader grouping and comparison of cell data.                                                                                                                                                                                                                                                                                                                                                                       |
| `obs["tissue_general_ontology_term_id"]`          | `string`  | (*Optional*) Ontology term identifier for the general tissue category, aiding in standardizing and grouping tissue types. For organoid or tissue samples, the Uber-anatomy ontology (`UBERON:`) is used. The term ids must be a child term of `UBERON:0001062` (anatomical entity). For cell cultures, the Cell Ontology (`CL:`) is used. The term ids cannot be `CL:0000255`, `CL:0000257` or `CL:0000548`.                                                                                  |
| `obs["batch"]`                                    | `string`  | (*Optional*) A batch identifier. This label is very context-dependent and may be a combination of the tissue, assay, donor, etc.                                                                                                                                                                                                                                                                                                                                                              |
| `obs["soma_joinid"]`                              | `integer` | (*Optional*) If the dataset was retrieved from CELLxGENE census, this is a unique identifier for the cell.                                                                                                                                                                                                                                                                                                                                                                                    |
| `obs["size_factors"]`                             | `double`  | (*Optional*) The size factors created by the normalisation method, if any.                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `var["feature_id"]`                               | `string`  | (*Optional*) Unique identifier for the feature, usually a ENSEMBL gene id.                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `var["feature_name"]`                             | `string`  | A human-readable name for the feature, usually a gene symbol.                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `var["soma_joinid"]`                              | `integer` | (*Optional*) If the dataset was retrieved from CELLxGENE census, this is a unique identifier for the feature.                                                                                                                                                                                                                                                                                                                                                                                 |
| `var["hvg"]`                                      | `boolean` | Whether or not the feature is considered to be a ‘highly variable gene’.                                                                                                                                                                                                                                                                                                                                                                                                                      |
| `var["hvg_score"]`                                | `integer` | A ranking of the features by hvg.                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| `obsm["X_pca"]`                                   | `double`  | The resulting PCA embedding.                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `obsp["knn_distances"]`                           | `double`  | K nearest neighbors distance matrix.                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| `obsp["knn_connectivities"]`                      | `double`  | K nearest neighbors connectivities matrix.                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `varm["pca_loadings"]`                            | `double`  | The PCA loadings matrix.                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| `layers["counts"]`                                | `integer` | Raw counts.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `layers["normalized"]`                            | `double`  | Normalised expression values.                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `uns["dataset_id"]`                               | `string`  | A unique identifier for the dataset. This is different from the `obs.dataset_id` field, which is the identifier for the dataset from which the cell data is derived.                                                                                                                                                                                                                                                                                                                          |
| `uns["dataset_name"]`                             | `string`  | A human-readable name for the dataset.                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| `uns["dataset_url"]`                              | `string`  | (*Optional*) Link to the original source of the dataset.                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| `uns["dataset_reference"]`                        | `string`  | (*Optional*) Bibtex reference of the paper in which the dataset was published.                                                                                                                                                                                                                                                                                                                                                                                                                |
| `uns["dataset_summary"]`                          | `string`  | Short description of the dataset.                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| `uns["dataset_description"]`                      | `string`  | Long description of the dataset.                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| `uns["dataset_organism"]`                         | `string`  | (*Optional*) The organism of the sample in the dataset.                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| `uns["normalization_id"]`                         | `string`  | Which normalization was used.                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `uns["pca_variance"]`                             | `double`  | The PCA variance objects.                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| `uns["knn"]`                                      | `object`  | Supplementary K nearest neighbors data.                                                                                                                                                                                                                                                                                                                                                                                                                                                       |

</div>

## Component type: Data processor

Path:
[`src/spatial_decomposition`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/spatial_decomposition)

A spatial decomposition dataset processor.

Arguments:

<div class="small">

| Name                      | Type   | Description                                                    |
|:--------------------------|:-------|:---------------------------------------------------------------|
| `--input`                 | `file` | A dataset processed by the common dataset processing pipeline. |
| `--output_single_cell`    | `file` | (*Output*) The single cell data file.                          |
| `--output_spatial_masked` | `file` | (*Output*) The masked spatial data file.                       |
| `--output_solution`       | `file` | (*Output*) The solution spatial data file.                     |

</div>

## File format: Single cell data

The single cell data file

Example file:
`resources_test/spatial_decomposition/pancreas/single_cell_ref.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obs: 'cell_type'
     layers: 'counts'
     uns: 'cell_type_names', 'dataset_id'

</div>

Slot description:

<div class="small">

| Slot                     | Type      | Description                                            |
|:-------------------------|:----------|:-------------------------------------------------------|
| `obs["cell_type"]`       | `integer` | Cell type label IDs.                                   |
| `layers["counts"]`       | `integer` | Raw counts.                                            |
| `uns["cell_type_names"]` | `string`  | Cell type names corresponding to values in `celltype`. |
| `uns["dataset_id"]`      | `string`  | A unique identifier for the dataset.                   |

</div>

## File format: Masked

The masked spatial data file

Example file:
`resources_test/spatial_decomposition/pancreas/spatial_masked.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obsm: 'coordinates'
     layers: 'counts'
     uns: 'cell_type_names', 'dataset_id'

</div>

Slot description:

<div class="small">

| Slot                     | Type      | Description                                                          |
|:-------------------------|:----------|:---------------------------------------------------------------------|
| `obsm["coordinates"]`    | `double`  | XY coordinates for each spot.                                        |
| `layers["counts"]`       | `integer` | Raw counts.                                                          |
| `uns["cell_type_names"]` | `string`  | Cell type names corresponding to columns of `proportions` in output. |
| `uns["dataset_id"]`      | `string`  | A unique identifier for the dataset.                                 |

</div>

## File format: Solution

The solution spatial data file

Example file:
`resources_test/spatial_decomposition/pancreas/solution.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obsm: 'coordinates', 'proportions_true'
     layers: 'counts'
     uns: 'cell_type_names', 'dataset_id'

</div>

Slot description:

<div class="small">

| Slot                       | Type      | Description                                                |
|:---------------------------|:----------|:-----------------------------------------------------------|
| `obsm["coordinates"]`      | `double`  | XY coordinates for each spot.                              |
| `obsm["proportions_true"]` | `double`  | True cell type proportions for each spot.                  |
| `layers["counts"]`         | `integer` | Raw counts.                                                |
| `uns["cell_type_names"]`   | `string`  | Cell type names corresponding to columns of `proportions`. |
| `uns["dataset_id"]`        | `string`  | A unique identifier for the dataset.                       |

</div>

## Component type: Control method

Path:
[`src/spatial_decomposition/control_methods`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/spatial_decomposition/control_methods)

Quality control methods for verifying the pipeline.

Arguments:

<div class="small">

| Name                     | Type   | Description                                         |
|:-------------------------|:-------|:----------------------------------------------------|
| `--input_single_cell`    | `file` | The single cell data file.                          |
| `--input_spatial_masked` | `file` | The masked spatial data file.                       |
| `--input_solution`       | `file` | The solution spatial data file.                     |
| `--output`               | `file` | (*Output*) Spatial data with estimated proportions. |

</div>

## Component type: Method

Path:
[`src/spatial_decomposition/methods`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/spatial_decomposition/methods)

A spatial composition method.

Arguments:

<div class="small">

| Name                  | Type   | Description                                         |
|:----------------------|:-------|:----------------------------------------------------|
| `--input_single_cell` | `file` | The single cell data file.                          |
| `--input_spatial`     | `file` | The masked spatial data file.                       |
| `--output`            | `file` | (*Output*) Spatial data with estimated proportions. |

</div>

## Component type: Metric

Path:
[`src/spatial_decomposition/metrics`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/spatial_decomposition/metrics)

A spatial decomposition metric.

Arguments:

<div class="small">

| Name               | Type   | Description                              |
|:-------------------|:-------|:-----------------------------------------|
| `--input_method`   | `file` | Spatial data with estimated proportions. |
| `--input_solution` | `file` | The solution spatial data file.          |
| `--output`         | `file` | (*Output*) Metric score file.            |

</div>

## File format: Output

Spatial data with estimated proportions.

Example file:
`resources_test/spatial_decomposition/pancreas/anndata_output.h5ad`

Description:

Spatial data file with estimated cell type proportions.

Format:

<div class="small">

    AnnData object
     obsm: 'coordinates', 'proportions_pred'
     layers: 'counts'
     uns: 'cell_type_names', 'dataset_id', 'method_id'

</div>

Slot description:

<div class="small">

| Slot                       | Type      | Description                                                |
|:---------------------------|:----------|:-----------------------------------------------------------|
| `obsm["coordinates"]`      | `double`  | XY coordinates for each spot.                              |
| `obsm["proportions_pred"]` | `double`  | Estimated cell type proportions for each spot.             |
| `layers["counts"]`         | `integer` | Raw counts.                                                |
| `uns["cell_type_names"]`   | `string`  | Cell type names corresponding to columns of `proportions`. |
| `uns["dataset_id"]`        | `string`  | A unique identifier for the dataset.                       |
| `uns["method_id"]`         | `string`  | A unique identifier for the method.                        |

</div>

## File format: Score

Metric score file.

Example file:
`resources_test/spatial_decomposition/pancreas/anndata_score.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     uns: 'dataset_id', 'method_id', 'metric_ids', 'metric_values'

</div>

Slot description:

<div class="small">

| Slot                   | Type     | Description                                                                                  |
|:-----------------------|:---------|:---------------------------------------------------------------------------------------------|
| `uns["dataset_id"]`    | `string` | A unique identifier for the dataset.                                                         |
| `uns["method_id"]`     | `string` | A unique identifier for the method.                                                          |
| `uns["metric_ids"]`    | `string` | One or more unique metric identifiers.                                                       |
| `uns["metric_values"]` | `double` | The metric values obtained for the given prediction. Must be of same length as ‘metric_ids’. |

</div>
