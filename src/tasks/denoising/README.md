# Denoising

Path:
[`src/tasks/denoising`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/src/tasks/denoising)

Single-cell RNA-Seq protocols only detect a fraction of the mRNA
molecules present in each cell. As a result, the measurements (UMI
counts) observed for each gene and each cell are associated with
generally high levels of technical noise ([Grün et al.,
2014](https://www.nature.com/articles/nmeth.2930)). Denoising describes
the task of estimating the true expression level of each gene in each
cell. In the single-cell literature, this task is also referred to as
*imputation*, a term which is typically used for missing data problems
in statistics. Similar to the use of the terms “dropout”, “missing
data”, and “technical zeros”, this terminology can create confusion
about the underlying measurement process ([Sarkar and Stephens,
2020](https://www.biorxiv.org/content/10.1101/2020.04.07.030007v2)).

A key challenge in evaluating denoising methods is the general lack of a
ground truth. A recent benchmark study ([Hou et al.,
2020](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02132-x))
relied on flow-sorted datasets, mixture control experiments ([Tian et
al., 2019](https://www.nature.com/articles/s41592-019-0425-8)), and
comparisons with bulk RNA-Seq data. Since each of these approaches
suffers from specific limitations, it is difficult to combine these
different approaches into a single quantitative measure of denoising
accuracy. Here, we instead rely on an approach termed molecular
cross-validation (MCV), which was specifically developed to quantify
denoising accuracy in the absence of a ground truth ([Batson et al.,
2019](https://www.biorxiv.org/content/10.1101/786269v1)). In MCV, the
observed molecules in a given scRNA-Seq dataset are first partitioned
between a *training* and a *test* dataset. Next, a denoising method is
applied to the training dataset. Finally, denoising accuracy is measured
by comparing the result to the test dataset. The authors show that both
in theory and in practice, the measured denoising accuracy is
representative of the accuracy that would be obtained on a ground truth
dataset.

``` mermaid
flowchart LR
  anndata_train(Training data)
  anndata_test(Test data)
  anndata_denoised(Denoised data)
  anndata_score(Score)
  anndata_common_dataset(NA)
  comp_control_method[/Control method/]
  comp_method[/Method/]
  comp_metric[/Metric/]
  comp_process_dataset[/Data processor/]
  anndata_train---comp_control_method
  anndata_test---comp_control_method
  anndata_train---comp_method
  anndata_test---comp_metric
  anndata_denoised---comp_metric
  anndata_common_dataset---comp_process_dataset
  comp_control_method-->anndata_denoised
  comp_method-->anndata_denoised
  comp_metric-->anndata_score
  comp_process_dataset-->anndata_train
  comp_process_dataset-->anndata_test
```

## File format: Common dataset

Example file: `resources_test/common/pancreas/dataset.h5ad`

A dataset processed by the common dataset processing pipeline. This
dataset contains both raw counts and normalized data matrices, as well
as a PCA embedding, HVG selection and a kNN graph.

Format:

<div class="small">

    AnnData object
     obs: 'celltype', 'batch', 'tissue', 'size_factors'
     var: 'hvg', 'hvg_score'
     obsm: 'X_pca'
     obsp: 'knn_distances', 'knn_connectivities'
     varm: 'pca_loadings'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_name', 'data_url', 'data_reference', 'dataset_summary', 'dataset_description', 'dataset_organism', 'pca_variance', 'knn'

</div>

Slot description:

<div class="small">

| Slot                         | Type      | Description                                                                    |
|:-----------------------------|:----------|:-------------------------------------------------------------------------------|
| `obs["celltype"]`            | `string`  | (*Optional*) Cell type information.                                            |
| `obs["batch"]`               | `string`  | (*Optional*) Batch information.                                                |
| `obs["tissue"]`              | `string`  | (*Optional*) Tissue information.                                               |
| `obs["size_factors"]`        | `double`  | (*Optional*) The size factors created by the normalisation method, if any.     |
| `var["hvg"]`                 | `boolean` | Whether or not the feature is considered to be a ‘highly variable gene’.       |
| `var["hvg_score"]`           | `integer` | A ranking of the features by hvg.                                              |
| `obsm["X_pca"]`              | `double`  | The resulting PCA embedding.                                                   |
| `obsp["knn_distances"]`      | `double`  | K nearest neighbors distance matrix.                                           |
| `obsp["knn_connectivities"]` | `double`  | K nearest neighbors connectivities matrix.                                     |
| `varm["pca_loadings"]`       | `double`  | The PCA loadings matrix.                                                       |
| `layers["counts"]`           | `integer` | Raw counts.                                                                    |
| `layers["normalized"]`       | `double`  | Normalised expression values.                                                  |
| `uns["dataset_id"]`          | `string`  | A unique identifier for the dataset.                                           |
| `uns["dataset_name"]`        | `string`  | Nicely formatted name.                                                         |
| `uns["data_url"]`            | `string`  | (*Optional*) Link to the original source of the dataset.                       |
| `uns["data_reference"]`      | `string`  | (*Optional*) Bibtex reference of the paper in which the dataset was published. |
| `uns["dataset_summary"]`     | `string`  | Short description of the dataset.                                              |
| `uns["dataset_description"]` | `string`  | Long description of the dataset.                                               |
| `uns["dataset_organism"]`    | `string`  | (*Optional*) The organism of the sample in the dataset.                        |
| `uns["pca_variance"]`        | `double`  | The PCA variance objects.                                                      |
| `uns["knn"]`                 | `object`  | Neighbors data.                                                                |

</div>

## Component type: Data processor

Path:
[`src/denoising`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/denoising)

Prepare a common dataset for the denoising task.

Arguments:

<div class="small">

| Name             | Type   | Description                                                                                                                                                                                                |
|:-----------------|:-------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--input`        | `file` | (*Optional*) A dataset processed by the common dataset processing pipeline. This dataset contains both raw counts and normalized data matrices, as well as a PCA embedding, HVG selection and a kNN graph. |
| `--output_train` | `file` | (*Optional, Output*) The training data.                                                                                                                                                                    |
| `--output_test`  | `file` | (*Optional, Output*) The test data.                                                                                                                                                                        |

</div>

## File format: train.h5ad

Example file: `resources_test/denoising/pancreas/train.h5ad`

The training data

Format:

<div class="small">

    AnnData object
     layers: 'counts'
     uns: 'dataset_id'

</div>

Slot description:

<div class="small">

| Slot                | Type      | Description                          |
|:--------------------|:----------|:-------------------------------------|
| `layers["counts"]`  | `integer` | Raw counts.                          |
| `uns["dataset_id"]` | `string`  | A unique identifier for the dataset. |

</div>

## File format: test.h5ad

Example file: `resources_test/denoising/pancreas/test.h5ad`

The test data

Format:

<div class="small">

    AnnData object
     layers: 'counts'
     uns: 'dataset_id'

</div>

Slot description:

<div class="small">

| Slot                | Type      | Description                          |
|:--------------------|:----------|:-------------------------------------|
| `layers["counts"]`  | `integer` | Raw counts.                          |
| `uns["dataset_id"]` | `string`  | A unique identifier for the dataset. |

</div>

## Component type: Control method

Path:
[`src/denoising/control_methods`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/denoising/control_methods)

This folder contains control components for the task. These components
have the same interface as the regular methods but also receive the
solution object as input. It serves as a starting point to test the
relative accuracy of new methods in the task, and also as a quality
control for the metrics defined in the task.

Arguments:

<div class="small">

| Name            | Type   | Description                             |
|:----------------|:-------|:----------------------------------------|
| `--input_train` | `file` | (*Optional*) The training data.         |
| `--input_test`  | `file` | (*Optional*) The test data.             |
| `--output`      | `file` | (*Optional, Output*) The denoised data. |

</div>

## Component type: Method

Path:
[`src/denoising/methods`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/denoising/methods)

A denoising method to remove noise (i.e. technical artifacts) from a
dataset.

Arguments:

<div class="small">

| Name            | Type   | Description                             |
|:----------------|:-------|:----------------------------------------|
| `--input_train` | `file` | (*Optional*) The training data.         |
| `--output`      | `file` | (*Optional, Output*) The denoised data. |

</div>

## Component type: Metric

Path:
[`src/denoising/metrics`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/denoising/metrics)

A metric for evaluating denoised datasets.

Arguments:

<div class="small">

| Name               | Type   | Description                             |
|:-------------------|:-------|:----------------------------------------|
| `--input_test`     | `file` | (*Optional*) The test data.             |
| `--input_denoised` | `file` | (*Optional*) The denoised data.         |
| `--output`         | `file` | (*Optional, Output*) Metric score file. |

</div>

## File format: magic.h5ad

Example file: `resources_test/denoising/pancreas/magic.h5ad`

The denoised data

Format:

<div class="small">

    AnnData object
     layers: 'counts', 'denoised'
     uns: 'dataset_id', 'method_id'

</div>

Slot description:

<div class="small">

| Slot                 | Type      | Description                          |
|:---------------------|:----------|:-------------------------------------|
| `layers["counts"]`   | `integer` | Raw counts.                          |
| `layers["denoised"]` | `integer` | denoised data.                       |
| `uns["dataset_id"]`  | `string`  | A unique identifier for the dataset. |
| `uns["method_id"]`   | `string`  | A unique identifier for the method.  |

</div>

## File format: magic_poisson.h5ad

Example file: `resources_test/denoising/pancreas/magic_poisson.h5ad`

Metric score file

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
