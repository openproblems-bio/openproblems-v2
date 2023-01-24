cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("batchelor", quietly = TRUE)
requireNamespace("SingleCellExperiment", quietly = TRUE)

## VIASH START
# path <- "resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter."
path <- "output/datasets/joint_embedding/openproblems_bmmc_multiome_phase1/openproblems_bmmc_multiome_phase1.censor_dataset.output_"
# path <- "output/public_datasets/joint_embedding/dyngen_citeseq_1/dyngen_citeseq_1.censor_dataset.output_"
par <- list(
  input_mod1 = paste0(path, "mod1.h5ad"),
  input_mod2 = paste0(path, "mod2.h5ad"),
  output = "output.h5ad",
  hvg_sel = 1000L
)
meta <- list(functionality_name = "foo")
## VIASH END

method_id <- meta$functionality_name

cat("Reading h5ad files\n")
input_mod1 <- anndata::read_h5ad(par$input_mod1)

rn <- rownames(input_mod1)
batch <- input_mod1$obs$batch
dataset_id <- input_mod1$uns[["dataset_id"]]
Xt_mod1 <- t(input_mod1$X)

# select hvg
if (!is.null(par$hvg_sel) && nrow(Xt_mod1) > par$hvg_sel) {
  sd_mod1 <- proxyC::rowSds(Xt_mod1)
  Xt_mod1 <- Xt_mod1[order(sd_mod1, decreasing = TRUE)[seq_len(par$hvg_sel)], ]
}

rm(input_mod1)
gc()

Xt_mod2 <- t(anndata::read_h5ad(par$input_mod2)$X)
if (!is.null(par$hvg_sel) && nrow(Xt_mod2) > par$hvg_sel) {
  sd_mod2 <- proxyC::rowSds(Xt_mod2)
  Xt_mod2 <- Xt_mod2[order(sd_mod2, decreasing = TRUE)[seq_len(par$hvg_sel)], ]
}

cat("Running fastMNN\n")
mnn_out <- batchelor::fastMNN(
  rbind(Xt_mod1, Xt_mod2),
  batch = batch
)
dr <- SingleCellExperiment::reducedDim(mnn_out, "corrected")

rownames(dr) <- rn
colnames(dr) <- paste0("comp_", seq_len(ncol(dr)))

out <- anndata::AnnData(
  X = dr,
  uns = list(
    dataset_id = dataset_id,
    method_id = meta$functionality_name
  ),
  obsm = list(X_emb = as(dr, "CsparseMatrix"))
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")
