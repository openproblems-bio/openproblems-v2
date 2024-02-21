# Cell-Cell Communication Source Target


Focuses on the prediction of interactions from steady-state, or
single-context, single-cell data between source cell types and target
cell types.

Path:
[`src/tasks/cell_cell_communication_source_target`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/tasks/cell_cell_communication_source_target)

## Motivation

The growing availability of single-cell data has sparked an increased
interest in the inference of cell-cell communication (CCC), with an
ever-growing number of computational tools developed for this purpose.

Different tools propose distinct preprocessing steps with diverse
scoring functions that are challenging to compare and evaluate.
Furthermore, each tool typically comes with its own set of prior
knowledge. To harmonize these, ([Dimitrov et al,
2022](https://openproblems.bio/documentation/reference/bibliography/#dimitrov2022comparison))
recently developed the ([LIANA+](https://github.com/saezlab/liana-py))
framework, which was used as a foundation for this task.

## Description

The challenges in evaluating the tools are further exacerbated by the
lack of a gold standard to benchmark the performance of CCC methods. In
an attempt to address this, Dimitrov et al use alternative data
modalities, including the spatial proximity of cell types and downstream
cytokine activities, to generate an inferred ground truth. However,
these modalities are only approximations of biological reality and come
with their own assumptions and limitations. In time, the inclusion of
more datasets with known ground truth interactions will become
available, from which the limitations and advantages of the different
CCC methods will be better understood.

This subtask evaluates methods in their ability to predict interactions
between spatially-adjacent source cell types and target cell types. This
subtask focuses on the prediction of interactions from steady-state, or
single-context, single-cell data.

## Authors & contributors

| name              | roles              |
|:------------------|:-------------------|
| Daniel Dimitrov   | maintainer, author |
| Scott Gigante     | author             |
| Daniel Strobl     | author             |
| Robrecht Cannoodt | contributor        |

## API

``` mermaid
flowchart LR
  file_common_dataset("Dataset")
  comp_process_dataset[/"Data processor"/]
  file_train("Training data")
  file_test("Test data")
  file_solution("Solution")
  comp_control_method[/"Control method"/]
  comp_method[/"Method"/]
  comp_metric[/"Metric"/]
  file_prediction("Prediction")
  file_score("Score")
  file_common_dataset---comp_process_dataset
  comp_process_dataset-->file_train
  comp_process_dataset-->file_test
  comp_process_dataset-->file_solution
  file_train---comp_control_method
  file_train---comp_method
  file_test---comp_control_method
  file_test---comp_method
  file_solution---comp_control_method
  file_solution---comp_metric
  comp_control_method-->file_prediction
  comp_method-->file_prediction
  comp_metric-->file_score
  file_prediction---comp_metric
```

## File format: Dataset

The dataset to pass to a method.

Example file:
`resources_test/cell_cell_communication_source_target/pancreas/dataset.h5ad`

Format:

<div class="small">

    AnnData object
     obs: 'label'
     var: 'hvg_score'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                         | Type      | Description                                                                          |
|:-----------------------------|:----------|:-------------------------------------------------------------------------------------|
| `obs["label"]`               | `string`  | Cell type annotation.                                                                |
| `var["hvg_score"]`           | `double`  | High variability gene score (normalized dispersion). The greater, the more variable. |
| `layers["counts"]`           | `integer` | Raw counts.                                                                          |
| `layers["normalized"]`       | `double`  | Normalized expression values.                                                        |
| `uns["dataset_id"]`          | `string`  | A unique identifier for the dataset.                                                 |
| `uns["dataset_name"]`        | `string`  | Nicely formatted name.                                                               |
| `uns["dataset_url"]`         | `string`  | (*Optional*) Link to the original source of the dataset.                             |
| `uns["dataset_reference"]`   | `string`  | (*Optional*) Bibtex reference of the paper in which the dataset was published.       |
| `uns["dataset_summary"]`     | `string`  | Short description of the dataset.                                                    |
| `uns["dataset_description"]` | `string`  | Long description of the dataset.                                                     |
| `uns["dataset_organism"]`    | `string`  | (*Optional*) The organism of the sample in the dataset.                              |
| `uns["normalization_id"]`    | `string`  | Which normalization was used.                                                        |

</div>

## Component type: Data processor

Path:
[`src/cell_cell_communication_source_target`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/cell_cell_communication_source_target)

A cell cell communication source_target dataset processor.

Arguments:

<div class="small">

| Name                | Type   | Description                                |
|:--------------------|:-------|:-------------------------------------------|
| `--input`           | `file` | The dataset to pass to a method.           |
| `--output_train`    | `file` | (*Output*) The training data.              |
| `--output_test`     | `file` | (*Output*) The test data.                  |
| `--output_solution` | `file` | (*Output*) The solution for the test data. |

</div>

## File format: Training data

The training data

Example file:
`resources_test/cell_cell_communication_source_target/pancreas/train.h5ad`

Format:

<div class="small">

    AnnData object
     obs: 'label'
     var: 'hvg_score'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                      | Type      | Description                          |
|:--------------------------|:----------|:-------------------------------------|
| `obs["label"]`            | `string`  | Cell type annotation.                |
| `var["hvg_score"]`        | `double`  | A ranking of the features by hvg.    |
| `layers["counts"]`        | `integer` | Raw counts.                          |
| `layers["normalized"]`    | `double`  | Normalized counts.                   |
| `uns["dataset_id"]`       | `string`  | A unique identifier for the dataset. |
| `uns["normalization_id"]` | `string`  | Which normalization was used.        |

</div>

## File format: Test data

The test data

Example file:
`resources_test/cell_cell_communication_source_target/pancreas/test.h5ad`

Format:

<div class="small">

    AnnData object
     obs: 'label'
     var: 'hvg_score'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                      | Type      | Description                          |
|:--------------------------|:----------|:-------------------------------------|
| `obs["label"]`            | `string`  | Cell type annotation.                |
| `var["hvg_score"]`        | `double`  | A ranking of the features by hvg.    |
| `layers["counts"]`        | `integer` | Raw counts.                          |
| `layers["normalized"]`    | `double`  | Normalized counts.                   |
| `uns["dataset_id"]`       | `string`  | A unique identifier for the dataset. |
| `uns["normalization_id"]` | `string`  | Which normalization was used.        |

</div>

## File format: Solution

The solution for the test data.

Example file:
`resources_test/cell_cell_communication_source_target/pancreas/solution.h5ad`

Format:

<div class="small">

    AnnData object
     obs: 'label'
     var: 'hvg_score'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                         | Type      | Description                                                                          |
|:-----------------------------|:----------|:-------------------------------------------------------------------------------------|
| `obs["label"]`               | `string`  | Cell type annotation.                                                                |
| `var["hvg_score"]`           | `double`  | High variability gene score (normalized dispersion). The greater, the more variable. |
| `layers["counts"]`           | `integer` | Raw counts.                                                                          |
| `layers["normalized"]`       | `double`  | Normalized expression values.                                                        |
| `uns["dataset_id"]`          | `string`  | A unique identifier for the dataset.                                                 |
| `uns["dataset_name"]`        | `string`  | Nicely formatted name.                                                               |
| `uns["dataset_url"]`         | `string`  | (*Optional*) Link to the original source of the dataset.                             |
| `uns["dataset_reference"]`   | `string`  | (*Optional*) Bibtex reference of the paper in which the dataset was published.       |
| `uns["dataset_summary"]`     | `string`  | Short description of the dataset.                                                    |
| `uns["dataset_description"]` | `string`  | Long description of the dataset.                                                     |
| `uns["dataset_organism"]`    | `string`  | (*Optional*) The organism of the sample in the dataset.                              |
| `uns["normalization_id"]`    | `string`  | Which normalization was used.                                                        |

</div>

## Component type: Control method

Path:
[`src/cell_cell_communication_source_target/control_methods`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/cell_cell_communication_source_target/control_methods)

Quality control methods for verifying the pipeline.

Arguments:

<div class="small">

| Name               | Type   | Description                     |
|:-------------------|:-------|:--------------------------------|
| `--input_train`    | `file` | The training data.              |
| `--input_test`     | `file` | The test data.                  |
| `--input_solution` | `file` | The solution for the test data. |
| `--output`         | `file` | (*Output*) The prediction file. |

</div>

## Component type: Method

Path:
[`src/cell_cell_communication_source_target/methods`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/cell_cell_communication_source_target/methods)

A Cell Cell Communication Source Target method.

Arguments:

<div class="small">

| Name            | Type   | Description                     |
|:----------------|:-------|:--------------------------------|
| `--input_train` | `file` | The training data.              |
| `--input_test`  | `file` | The test data.                  |
| `--output`      | `file` | (*Output*) The prediction file. |

</div>

## Component type: Metric

Path:
[`src/cell_cell_communication_source_target/metrics`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/cell_cell_communication_source_target/metrics)

A CCC source target metric.

Arguments:

<div class="small">

| Name                 | Type   | Description                     |
|:---------------------|:-------|:--------------------------------|
| `--input_solution`   | `file` | The solution for the test data. |
| `--input_prediction` | `file` | The prediction file.            |
| `--output`           | `file` | (*Output*) Metric score file.   |

</div>

## File format: Prediction

The prediction file

Example file:
`resources_test/cell_cell_communication_source_target/pancreas/prediction.h5ad`

Format:

<div class="small">

    AnnData object
     obs: 'ccc_pred'
     uns: 'dataset_id', 'normalization_id', 'method_id'

</div>

Slot description:

<div class="small">

| Slot                      | Type     | Description                                |
|:--------------------------|:---------|:-------------------------------------------|
| `obs["ccc_pred"]`         | `string` | Predicted interactions for the test cells. |
| `uns["dataset_id"]`       | `string` | A unique identifier for the dataset.       |
| `uns["normalization_id"]` | `string` | Which normalization was used.              |
| `uns["method_id"]`        | `string` | A unique identifier for the method.        |

</div>

## File format: Score

Metric score file

Example file:
`resources_test/cell_cell_communication_source_target/pancreas/score.h5ad`

Format:

<div class="small">

    AnnData object
     uns: 'dataset_id', 'normalization_id', 'method_id', 'metric_ids', 'metric_values'

</div>

Slot description:

<div class="small">

| Slot                      | Type     | Description                                                                                  |
|:--------------------------|:---------|:---------------------------------------------------------------------------------------------|
| `uns["dataset_id"]`       | `string` | A unique identifier for the dataset.                                                         |
| `uns["normalization_id"]` | `string` | Which normalization was used.                                                                |
| `uns["method_id"]`        | `string` | A unique identifier for the method.                                                          |
| `uns["metric_ids"]`       | `string` | One or more unique metric identifiers.                                                       |
| `uns["metric_values"]`    | `double` | The metric values obtained for the given prediction. Must be of same length as ‘metric_ids’. |

</div>

