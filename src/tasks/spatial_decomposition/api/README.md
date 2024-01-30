# Spatial Decomposition/Deconvolution

Calling cell-type compositions for spot-based spatial transcriptomics data

Path:
[`src/tasks/spatial_decomposition`]

## Motivation

Spatial decomposition (also often referred to as Spatial deconvolution) is
applicable to spatial transcriptomics data where the transcription profile of
each capture location (spot, voxel, bead, etc.) do not share a bijective
relationship with the cells in the tissue, i.e., multiple cells may contribute
to the same capture location. The task of spatial decomposition then refers to
estimating the composition of cell types/states that are present at each capture
location. The cell type/states estimates are presented as proportion values,
representing the proportion of the cells at each capture location that belong to
a given cell type.

We distinguish between _reference-based_ decomposition and _de novo_
decomposition, where the former leverage external data (e.g., scRNA-seq or
scNuc-seq) to guide the inference process, while the latter only work with the
spatial data. We require that all datasets have an associated reference single
cell data set, but methods are free to ignore this information.

## Description



## Authors & contributors



## API

``` mermaid
flowchart LR
  file_common_dataset("Common dataset")
  comp_process_dataset[/"Data processor"/]
  file_single_cell("Single cell reference")
  file_spatial_masked("Spatial data")
  file_solution("Solution")
  comp_control_method[/"Control method"/]
  comp_method[/"Method"/]
  comp_metric[/"Metric"/]
  file_output("Spatially Deconvoluted")
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

## Component type: Data processor

## File format: Dataset

## File format: Test data

## Component type: Control method

## Component type: Method

## Component type: Metric

## File format: Embedding

## File format: Score
