## VIASH START
par <- list(
  input_mod1 = "resources_test/common/multimodal/dataset_mod1.h5ad",
  input_mod2 = "resources_test/common/multimodal/dataset_mod2.h5ad",
  output = "output.mnn.h5ad"
)
## VIASH END

cat("Loading dependencies\n")
library(anndata, warn.conflicts = FALSE)
library(Matrix, warn.conflicts = FALSE)
requireNamespace("sparsesvd", quietly = TRUE)
requireNamespace("batchelor", quietly = TRUE)

cat("Reading input h5ad file\n")
adata_mod1 <- read_h5ad(par$input_mod1)
adata_mod2 <- read_h5ad(par$input_mod2)


cat("Running SVD\n")
mode1_svd <- adata_mod1$obsm[["svd"]]
mode1_svd_uv <- mode1_svd$u %*% diag(mode1_svd$d)
mode2_svd <- sparsesvd::sparsesvd(mode2, rank = n_svd)
mode2_svd_uv <- mode2_svd$u %*% diag(mode2_svd$d)

cat("Running MNN\n")
sce_mnn <- batchelor::fastMNN(
  t(mode1_svd_uv),
  t(mode2_svd_uv)
)

cat("Storing output\n")
combined_recons <- t(SummarizedExperiment::assay(sce_mnn, "reconstructed"))
mode1_recons <- combined_recons[seq_len(nrow(mode1_svd_uv)), , drop = FALSE]
mode2_recons <- combined_recons[-seq_len(nrow(mode1_svd_uv)), , drop = FALSE]

adata$obsm[["aligned"]] <- as.matrix(mode1_recons)
adata$obsm[["mode2_aligned"]] <- as.matrix(mode2_recons)

cat("Writing to file\n")
adata$uns["method_id"] = "mnn"
zzz <- adata$write_h5ad(par$output, compression = "gzip")
