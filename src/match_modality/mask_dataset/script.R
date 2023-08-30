cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(assertthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)

## VIASH START
# input_path <- "resources_test/common/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter."
# input_path <- "output/datasets/common/openproblems_bmmc_multiome_phase1/openproblems_bmmc_multiome_phase1.manual_formatting."
# input_path <- "output/datasets/common/openproblems_bmmc_cite_phase1/openproblems_bmmc_cite_phase1.manual_formatting."
input_path <- "resources_test/common/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter."
output_path <- "output/multiome"
# output_path <- "output/datasets/match_modality/openproblems_bmmc_multiome_phase1/openproblems_bmmc_multiome_phase1.censor_dataset."
# output_path <- "output/datasets/match_modality/openproblems_bmmc_multiome_iid/openproblems_bmmc_multiome_iid.censor_dataset."
# dir.create(dirname(output_path), recursive = TRUE)

par <- list(
  input_mod1 = paste0(input_path, "output_rna.h5ad"),
  input_mod2 = paste0(input_path, "output_mod2.h5ad"),
  output_train_mod1 = paste0(output_path, "output_train_mod1.h5ad"),
  output_train_mod2 = paste0(output_path, "output_train_mod2.h5ad"),
  output_train_sol = paste0(output_path, "output_train_sol.h5ad"),
  output_test_mod1 = paste0(output_path, "output_test_mod1.h5ad"),
  output_test_mod2 = paste0(output_path, "output_test_mod2.h5ad"),
  output_test_sol = paste0(output_path, "output_test_sol.h5ad"),
  seed = 1L,
  knn = 10L
)
## VIASH END

set.seed(par$seed)

subset_mats <- function(li, obs_filt, anonymize = FALSE) {
  out <- list()
  for (n in names(li)) {
    mat <- li[[n]][obs_filt, , drop = FALSE]
    if (anonymize) {
      rownames(mat) <- paste0("cell_", seq_len(nrow(mat)))
    }
    out[[n]] <- mat
  }
  out
}


cat("Reading input data\n")
input_mod1 <- anndata::read_h5ad(par$input_mod1)
input_mod2 <- anndata::read_h5ad(par$input_mod2)
ad1_mod <- unique(input_mod1$var[["feature_types"]])
ad2_mod <- unique(input_mod2$var[["feature_types"]])
new_dataset_id <- paste0(input_mod1$uns[["dataset_id"]], "_MM_", tolower(ad1_mod), "2", tolower(ad2_mod))
ad1_uns <- list(dataset_id = new_dataset_id, organism = "human")
ad2_uns <- list(dataset_id = new_dataset_id, organism = "human")
ad1_obsm <- list()
ad2_obsm <- list()

if (ad1_mod == "ATAC") {
  ad1_uns$gene_activity_var_names <- input_mod1$uns$gene_activity_var_names
  ad1_obsm$gene_activity <- as(input_mod1$obsm$gene_activity, "CsparseMatrix")
}
if (ad2_mod == "ATAC") {
  ad2_uns$gene_activity_var_names <- input_mod2$uns$gene_activity_var_names
  ad2_obsm$gene_activity <- as(input_mod2$obsm$gene_activity, "CsparseMatrix")
}

cat("Shuffle train cells\n")
train_ix <- which(input_mod1$obs$is_train) %>% sort
train_mod2_ix <- sample.int(length(train_ix))

cat("Shuffle test cells\n")
test_ix <- which(!input_mod1$obs$is_train) %>% sort
test_mod2_ix <- sample.int(length(test_ix))

is_categorical <- function(x) is.character(x) || is.factor(x)
# relevel <- function(x) factor(as.character(x))
relevel <- function(x) as.character(x)

cat("Creating train objects\n")
mod1_var <- input_mod1$var %>% select(one_of("gene_ids", "feature_types"))
mod2_var <- input_mod2$var %>% select(one_of("gene_ids", "feature_types"))
train_obs1 <- input_mod1$obs[train_ix, , drop = FALSE] %>%
  select(one_of("batch", "size_factors")) %>%
  mutate_if(is_categorical, relevel)
train_obs2 <- input_mod2$obs[train_ix[train_mod2_ix], , drop = FALSE] %>%
  select(one_of("batch", "size_factors")) %>%
  mutate_if(is_categorical, relevel)
rownames(train_obs2) <- NULL
if (ncol(train_obs2) == 0) train_obs2 <- NULL
assert_that("size_factors" %in% colnames(train_obs1) != "size_factors" %in% colnames(train_obs2))
assert_that(all(train_obs1$batch == train_obs2$batch[order(train_mod2_ix)]))

output_train_mod1 <- anndata::AnnData(
  X = input_mod1$X[train_ix, , drop = FALSE],
  layers = subset_mats(input_mod1$layers, train_ix),
  obsm = subset_mats(ad1_obsm, train_ix),
  obs = train_obs1,
  var = mod1_var,
  uns = ad1_uns
)
output_train_mod2 <- anndata::AnnData(
  X = input_mod2$X[train_ix[train_mod2_ix], , drop = FALSE] %>%
    magrittr::set_rownames(., paste0("cell_", seq_len(nrow(.)))),
  layers = subset_mats(input_mod2$layers, train_ix[train_mod2_ix], anonymize = TRUE),
  obsm = subset_mats(ad2_obsm, train_ix[train_mod2_ix], anonymize = TRUE),
  obs = train_obs2,
  var = mod2_var,
  uns = ad2_uns
)
assert_that(all(output_train_mod1$obs$batch == output_train_mod2$obs$batch[order(train_mod2_ix)]))

cat("Create test objects\n")
test_obs1 <- input_mod1$obs[test_ix, , drop = FALSE] %>%
  select(one_of("batch", "size_factors")) %>%
  mutate_if(is_categorical, relevel)
test_obs2 <- input_mod2$obs[test_ix[test_mod2_ix], , drop = FALSE] %>%
  select(one_of("batch", "size_factors")) %>%
  mutate_if(is_categorical, relevel)
rownames(test_obs2) <- NULL
if (ncol(test_obs2) == 0) test_obs2 <- NULL
assert_that("size_factors" %in% colnames(train_obs1) != "size_factors" %in% colnames(train_obs2))
assert_that(all(test_obs1$batch == test_obs2$batch[order(test_mod2_ix)]))

output_test_mod1 <- anndata::AnnData(
  X = input_mod1$X[test_ix, , drop = FALSE],
  layers = subset_mats(input_mod1$layers, test_ix),
  obsm = subset_mats(ad1_obsm, test_ix),
  obs = test_obs1,
  var = mod1_var,
  uns = ad1_uns
)
output_test_mod2 <- anndata::AnnData(
  X = input_mod2$X[test_ix[test_mod2_ix], , drop = FALSE] %>%
    magrittr::set_rownames(., paste0("cell_", seq_len(nrow(.)))),
  layers = subset_mats(input_mod2$layers, test_ix[test_mod2_ix], anonymize = TRUE),
  obsm = subset_mats(ad2_obsm, test_ix[test_mod2_ix], anonymize = TRUE),
  obs = test_obs2,
  var = mod2_var,
  uns = ad2_uns
)
assert_that(all(output_test_mod1$obs$batch == output_test_mod2$obs$batch[order(test_mod2_ix)]))

cat("Create solution objects\n")

train_sol_mat <- Matrix::sparseMatrix(
  i = seq_along(train_mod2_ix),
  j = order(train_mod2_ix),
  x = rep(1, length(train_mod2_ix))
)
output_train_sol <- anndata::AnnData(
  X = train_sol_mat,
  obs = input_mod1$obs[train_ix, , drop = FALSE] %>% select(one_of(c("batch"))) %>% mutate_if(is_categorical, relevel),
  uns = list(dataset_id = new_dataset_id, pairing_ix = train_mod2_ix - 1)
)

test_sol_mat <- Matrix::sparseMatrix(
  i = seq_along(test_mod2_ix),
  j = order(test_mod2_ix),
  x = rep(1, length(test_mod2_ix))
)
output_test_sol <- anndata::AnnData(
  X = test_sol_mat,
  obs = input_mod1$obs[test_ix, , drop = FALSE] %>% select(one_of(c("batch"))) %>% mutate_if(is_categorical, relevel),
  uns = list(dataset_id = new_dataset_id, pairing_ix = test_mod2_ix - 1)
)

# checks
# mean(rowSums(train_solknn > 0))
# mean(rowSums(test_solknn > 0))
# sum(train_solknn * train_sol_mat) == nrow(train_sol_mat)
# sum(test_solknn * test_sol_mat) == nrow(test_sol_mat)

cat("Saving output files as h5ad\n")
zzz <- output_train_mod1$write_h5ad(par$output_train_mod1, compression = "gzip")
zzz <- output_train_mod2$write_h5ad(par$output_train_mod2, compression = "gzip")
zzz <- output_train_sol$write_h5ad(par$output_train_sol, compression = "gzip")
zzz <- output_test_mod1$write_h5ad(par$output_test_mod1, compression = "gzip")
zzz <- output_test_mod2$write_h5ad(par$output_test_mod2, compression = "gzip")
zzz <- output_test_sol$write_h5ad(par$output_test_sol, compression = "gzip")
