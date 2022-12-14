cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(assertthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)

## VIASH START
input_path <- "output/datasets_2021-11-08/common/openproblems_bmmc_multiome_phase1v2/openproblems_bmmc_multiome_phase1v2.manual_formatting."
output_path <- ""

par <- list(
  input_mod1 = paste0(input_path, "output_rna.h5ad"),
  input_mod2 = paste0(input_path, "output_mod2.h5ad"),
  output_mod1 = paste0(output_path, "output_mod1.h5ad"),
  output_mod2 = paste0(output_path, "output_mod2.h5ad"),
  output_solution = paste0(output_path, "solution.h5ad"),
  train_only = TRUE
)
## VIASH END

cat("Reading mod1 data\n")
input_mod1 <- anndata::read_h5ad(par$input_mod1)
ad1_mod <- unique(input_mod1$var[["feature_types"]])
new_dataset_id <- paste0(input_mod1$uns[["dataset_id"]], "_JE")
ad1_uns <- list(dataset_id = new_dataset_id, organism = "human")
ad2_uns <- list(dataset_id = new_dataset_id, organism = "human")

cat("Creating mod1 object\n")
out_mod1 <- anndata::AnnData(
  X = input_mod1$X,
  layers = list(counts = input_mod1$layers[["counts"]]),
  var = input_mod1$var %>% select(one_of("gene_ids"), feature_types),
  obs = input_mod1$obs %>% select(one_of("batch", "size_factors")),
  uns = ad1_uns
)

cat("Create solution object\n")
out_solution <- anndata::AnnData(
  X = input_mod1$X,
  var = input_mod1$var %>% select(one_of("gene_ids"), feature_types),
  obs = input_mod1$obs %>% select(
    one_of("batch", "cell_type", "pseudotime_order_GEX", "pseudotime_order_ATAC", "pseudotime_order_ADT", "S_score", "G2M_score")
  ),
  uns = ad1_uns
)

is_train <- input_mod1$obs$is_train

if (par$train_only) {
  cat("Filtering out test cells\n", sep = "")
  out_mod1 <- out_mod1[is_train, ] #$copy()
  out_solution <- out_solution[is_train, ]# $copy()
}

rm(input_mod1)
gc()

cat("Reading mod2 data\n")
input_mod2 <- anndata::read_h5ad(par$input_mod2)
ad2_mod <- unique(input_mod2$var[["feature_types"]])
ad2_obsm <- list()

if (ad2_mod == "ATAC") {
  ad2_uns$gene_activity_var_names <- input_mod2$uns$gene_activity_var_names
  ad2_obsm$gene_activity <- as(input_mod2$obsm$gene_activity, "CsparseMatrix")
}

cat("Creating mod2 object\n")
out_mod2 <- anndata::AnnData(
  X = input_mod2$X,
  layers = list(counts = input_mod2$layers[["counts"]]),
  var = input_mod2$var %>% select(one_of("gene_ids"), feature_types),
  obs = input_mod2$obs %>% select(one_of("batch")),
  obsm = ad2_obsm,
  uns = ad2_uns
)
rm(input_mod2)
gc()

if (par$train_only) {
  cat("Filtering out test cells\n", sep = "")
  out_mod2 <- out_mod2[is_train, ] #$copy()
}

cat("Saving output files as h5ad\n")
cat("output_mod1:")
print(out_mod1)
zzz <- out_mod1$write_h5ad(par$output_mod1, compression = "gzip")

cat("output_mod2:")
print(out_mod2)
zzz <- out_mod2$write_h5ad(par$output_mod2, compression = "gzip")

cat("output_solution:")
print(out_solution)
zzz <- out_solution$write_h5ad(par$output_solution, compression = "gzip")
