cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("NewWave", quietly = TRUE)
requireNamespace("SingleCellExperiment", quietly = TRUE)

## VIASH START
# path <- "resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter."
path <- "output/datasets/joint_embedding/openproblems_bmmc_multiome_phase1/openproblems_bmmc_multiome_phase1.censor_dataset.output_"
# path <- "output/public_datasets/joint_embedding/dyngen_citeseq_1/dyngen_citeseq_1.censor_dataset.output_"
par <- list(
  input_mod1 = paste0(path, "mod1.h5ad"),
  input_mod2 = paste0(path, "mod2.h5ad"),
  output = "output.h5ad",
  maxiter = 2L,
  k = 3L,
  hvg_sel = 1000
)
meta <- list(functionality_name = "foo")
## VIASH END

method_id <- meta$functionality_name

cat("Reading mod1 h5ad\n")
input_mod1 <- anndata::read_h5ad(par$input_mod1)

rn <- rownames(input_mod1)
batch <- input_mod1$obs$batch
dataset_id <- input_mod1$uns[["dataset_id"]]

sd1 <- proxyC::colSds(input_mod1$X)
fil1 <-
  if (!is.null(par$hvg_sel) && ncol(input_mod1) > par$hvg_sel) {
    head(order(sd1, decreasing = TRUE), par$hvg_sel)
  } else {
    which(sd1 > 0)
  }
data1 <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = t(input_mod1$layers[["counts"]][, fil1])),
  colData = data.frame(batch = factor(batch))
)
rm(input_mod1)
gc()

cat("Running NewWave on mod1\n")
res1 <- NewWave::newWave(
  data1,
  X = "~batch",
  verbose = TRUE,
  K = par$k,
  maxiter_optimize = par$maxiter,
  n_gene_par = min(300, nrow(data1)),
  n_cell_par = min(300, ncol(data1)),
  commondispersion = FALSE
)
rm(data1)

dr_x1 <- SingleCellExperiment::reducedDim(res1)

cat("Reading mod2 anndata\n")
input_mod2 <- anndata::read_h5ad(par$input_mod2)
sd2 <- proxyC::colSds(input_mod2$X)
fil2 <-
  if (!is.null(par$hvg_sel) && ncol(input_mod2) > par$hvg_sel) {
    head(order(sd2, decreasing = TRUE), par$hvg_sel)
  } else {
    which(sd2 > 0)
  }
data2 <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = t(input_mod2$layers[["counts"]][, fil2])),
  colData = data.frame(batch = factor(batch))
)
rm(input_mod2)
gc()

cat("Running NewWave on mod2\n")
res2 <- NewWave::newWave(
  data2,
  X = "~batch",
  verbose = TRUE,
  K = par$k,
  maxiter_optimize = par$maxiter,
  n_gene_par = min(300, nrow(data2)),
  n_cell_par = min(300, ncol(data2)),
  commondispersion = FALSE
)
dr_x2 <- SingleCellExperiment::reducedDim(res2)
rm(data2)

cat("Spline separate DRs\n")
dr <- do.call(cbind, lapply(seq_len(ncol(dr_x1)), function(i) {
  cbind(dr_x1[, i], dr_x2[, i])
}))

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
