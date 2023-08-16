cat("Loading dependencies\n")
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE)
requireNamespace("batchelor", quietly = TRUE)
library(SingleCellExperiment, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input = 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
  output = 'output.h5ad',
  hvg = FALSE
)
## VIASH END

cat("Read input\n")
adata <- anndata::read_h5ad(par$input)

# don't subset when return_type is not "feature"
if ("hvg" %in% names(par) && par$hvg) {
  cat("Select HVGs\n")
  adata <- adata[, adata$var[["hvg"]]]
}

cat("Run mnn\n")
out <- batchelor::fastMNN(
  t(adata$layers[["normalized"]]),
  batch = adata$obs[["batch"]]
)

cat("Store output\n")
# reusing the same script for mnn_correct and mnn_correct_feature
return_type <- gsub("mnn_correct_", "", meta[["functionality_name"]])

if (return_type == "feature") {
  layer <- SummarizedExperiment::assay(out, "reconstructed")
  adata$layers[["corrected_counts"]] <- as(t(layer), "sparseMatrix")
} else if (return_type == "embedding") {
  obsm <- SingleCellExperiment::reducedDim(out, "corrected")
  adata$obsm[["X_emb"]] <- obsm
}

cat("Store outputs\n")
adata$uns[["method_id"]] <- meta$functionality_name
adata$write_h5ad(par$output, compression = "gzip")
