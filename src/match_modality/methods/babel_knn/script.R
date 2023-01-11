cat(">> Loading dependencies\n")

requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
library(testthat, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("DropletUtils", quietly = TRUE)

options(tidyverse.quiet = TRUE)
library(tidyverse)

babel_location <- "/babel/bin/"
conda_bin <- "/opt/conda/bin/conda"

## VIASH START
path <- "output/datasets/match_modality/openproblems_bmmc_multiome_phase1/openproblems_bmmc_multiome_phase1.censor_dataset.output_"
par <- list(
  input_train_mod1 = paste0(path, "train_mod1.h5ad"),
  input_train_mod2 = paste0(path, "train_mod2.h5ad"),
  input_train_sol = paste0(path, "train_sol.h5ad"),
  input_test_mod1 = paste0(path, "test_mod1.h5ad"),
  input_test_mod2 = paste0(path, "test_mod2.h5ad"),
  output = "output.h5ad",
  n_dims = 10,
  n_neighs = 10
)
conda_bin <- "conda"
babel_location <- "../babel/bin/"
## VIASH END

input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2, backed = TRUE)
if (input_train_mod2$var$feature_types[[1]] != "ATAC") {
  cat("Error: babel only runs on GEX to ATAC datasets\n")
  quit(save = "no", status = 42)
}

cat("Reading h5ad files\n")
input_train_mod1 <- anndata::read_h5ad(par$input_train_mod1)
input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2)
input_train_sol <- anndata::read_h5ad(par$input_train_sol)
input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1)
input_test_mod2 <- anndata::read_h5ad(par$input_test_mod2)

cat(">> Reading h5ad files\n")
if (is.null(input_train_mod1$var$gene_ids)) input_train_mod1$var$gene_ids <- colnames(input_train_mod1)
if (is.null(input_train_mod2$var$gene_ids)) input_train_mod2$var$gene_ids <- colnames(input_train_mod2)
if (is.null(input_test_mod1$var$gene_ids)) input_test_mod1$var$gene_ids <- colnames(input_test_mod1)
if (is.null(input_test_mod2$var$gene_ids)) input_test_mod2$var$gene_ids <- colnames(input_test_mod2)

mod1 <- as.character(unique(input_train_mod1$var$feature_types))
mod2 <- as.character(unique(input_train_mod2$var$feature_types))

# multiome_matrix for export to Babel's input format
multiome_matrix <- cbind(input_train_mod1$X, input_train_mod2$X)

# generate multiome anndata object for training
ad_babel <- anndata::AnnData(
  X = multiome_matrix,
  var = bind_rows(input_train_mod1$var, input_train_mod2$var)
)

# setting up babel dirs
tmpdir <- tempfile(pattern = "babel_temp", fileext = "/")
cat(">> Setting up directories for babel at ", tmpdir, "\n", sep = "")
dir.create(tmpdir)
on.exit(unlink(tmpdir, recursive = TRUE))

dir_data <- paste0(tmpdir, "data/")     # location of input files
dir.create(dir_data)
dir_model <- paste0(tmpdir, "model/")   # location of babel model
dir_pred <- paste0(tmpdir, "pred/")     # location of predictions

feature_type_map <- c(
  "GEX" = "Gene Expression",
  "ADT" = "Peaks", # try to make it run on ADT data as well
  "ATAC" = "Peaks"
)

cat(">> Writing train dataset as 10x-CellRanger H5 format\n")
DropletUtils::write10xCounts(
  paste0(dir_data, "train_input.h5"),
  t(ad_babel$X),
  gene.id = ad_babel$var$gene_ids,
  gene.symbol = colnames(ad_babel),
  barcodes = rownames(ad_babel),
  type = "HDF5",
  version = "3",
  genome = "GRCh38",
  gene.type = unname(feature_type_map[ad_babel$var$feature_types]),
  overwrite = TRUE
)

cat(">> Writing test dataset as 10x-CellRanger H5 format\n")
DropletUtils::write10xCounts(
  paste0(dir_data, "test_input.h5"),
  t(input_test_mod1$X),
  gene.id = input_test_mod1$var$gene_ids,
  gene.symbol = colnames(input_test_mod1),
  barcodes = rownames(input_test_mod1),
  type = "HDF5",
  version = "3",
  genome = "GRCh38",
  gene.type = unname(feature_type_map[input_test_mod1$var$feature_types]),
  overwrite = TRUE
)


cat(">> Babel: train model\n")
babel_train_cmd <- paste0(
  conda_bin, " run -n babel ",
  "python ", babel_location, "train_model.py ",
  "--data ", dir_data, "train_input.h5 ",
  "--outdir ", dir_model, " ",
  "--nofilter"
)
# stringent filtering causes babel to sometimes fail
# reason: https://github.com/wukevin/babel/blob/main/babel/sc_data_loaders.py#L168-L190

out1 <- system(babel_train_cmd)

# check whether training succeeded
expect_equal(out1, 0, info = paste0("Model training failed with exit code ", out1))

cat(">> Babel: predict from model\n")
babel_pred_cmd <- paste0(
  conda_bin, " run -n babel ",
  "python ", babel_location, "predict_model.py ",
  "--checkpoint ", dir_model, " ",
  "--data ", dir_data, "test_input.h5 ",
  "--outdir ", dir_pred, " ",
  "--nofilter"
)
out2 <- system(babel_pred_cmd)

# check whether prediction succeeded
expect_equal(out2, 0, info = paste0("Prediction failed with exit code ", out1))

cat(">> Read predictions\n")
pred <- anndata::read_h5ad(paste0(dir_pred, "/rna_atac_adata.h5ad"))

#######################################
#####  KNN

# Only some features are present in Babel's output
input_test_mod2_filter <- input_test_mod2[, row.names(input_test_mod2$var) %in% row.names(pred$var)]

# Dimensional reduction of both predicted and test profiles
pred_profiles <- pred[, row.names(input_test_mod2_filter$var)]$X

cat("Performing DR on test values\n")
dr <- lmds::lmds(
  rbind(pred_profiles, input_test_mod2_filter$X),
  ndim = par$n_dims,
  distance_method = par$distance_method
)

train_ix <- seq_len(nrow(pred_profiles))
dr_preds <- dr[train_ix, , drop = FALSE]
dr_test_mod2 <- dr[-train_ix, , drop = FAPSE]


cat("Performing KNN between test mod2 DR and predicted test mod2\n")
knn_out <- FNN::get.knnx(
  dr_preds,
  dr_test_mod2,
  k = min(1000, nrow(dr_preds))
)

cat("Creating output data structures\n")
df <- tibble(
  i = as.vector(row(knn_out$nn.index)),
  j = as.vector(knn_out$nn.index),
  x = max(knn_out$nn.dist) * 2 - as.vector(knn_out$nn.dist)
)
knn_mat <- Matrix::sparseMatrix(
  i = df$i,
  j = df$j,
  x = df$x,
  dims = list(nrow(input_test_mod1), nrow(input_test_mod2))
)

# normalise to make rows sum to 1
rs <- Matrix::rowSums(knn_mat)
knn_mat@x <- knn_mat@x / rs[knn_mat@i + 1]

cat("Creating output anndata\n")
out <- anndata::AnnData(
  X = as(knn_mat, "CsparseMatrix"),
  uns = list(
    dataset_id = input_train_mod1$uns[["dataset_id"]],
    method_id = "baseline_babel_knn"
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")