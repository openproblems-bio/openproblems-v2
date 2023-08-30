cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
path <- "output/datasets/match_modality/openproblems_bmmc_multiome_phase1_rna/openproblems_bmmc_multiome_phase1_rna.censor_dataset.output_"
path <- "output/datasets/match_modality/openproblems_bmmc_cite_phase1_rna/openproblems_bmmc_cite_phase1_rna.censor_dataset.output_"
par <- list(
  input_train_mod1 = paste0(path, "train_mod1.h5ad"),
  input_train_mod2 = paste0(path, "train_mod2.h5ad"),
  input_train_sol = paste0(path, "train_sol.h5ad"),
  input_test_mod1 = paste0(path, "test_mod1.h5ad"),
  input_test_mod2 = paste0(path, "test_mod2.h5ad"),
  output = "output.h5ad",
  n_pop = 300L
)
meta <- list(functionality_name = "foo")
## VIASH END

n_cores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)

method_id <- meta$functionality_name

cat("Read train sol\n")
input_train_sol <- anndata::read_h5ad(par$input_train_sol)

cat("Reading mod1 h5ad files\n")
input_train_mod1 <- anndata::read_h5ad(par$input_train_mod1)
input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1)

cat("Running LMDS on input data\n")
# merge input matrices
mod1_X <- rbind(input_train_mod1$X, input_test_mod1$X)
rm(input_train_mod1, input_test_mod1)
gc()

# perform DR
dr_x1 <- lmds::lmds(mod1_X, ndim = 10, distance_method = "pearson")
rm(mod1_X)
gc()

# split input matrices
dr_x1_train <- dr_x1[seq_len(nrow(input_train_sol)), , drop = FALSE]
dr_x1_test <- dr_x1[-seq_len(nrow(input_train_sol)), , drop = FALSE]

cat("Reading mod1 h5ad files\n")
input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2)
input_test_mod2 <- anndata::read_h5ad(par$input_test_mod2)

cat("Running LMDS on input data\n")
# merge input matrices
match_train <- input_train_sol$uns$pairing_ix + 1
mod2_X <- rbind(input_train_mod2$X[order(match_train), , drop = FALSE], input_test_mod2$X)
rm(input_train_mod2, input_test_mod2)
gc()

# perform DR
dr_x2 <- lmds::lmds(mod2_X, ndim = 3, distance_method = "pearson")
rm(mod2_X)
gc()

# split input matrices
dr_x2_train <- dr_x2[seq_len(nrow(input_train_sol)), , drop = FALSE]
dr_x2_test <- dr_x2[-seq_len(nrow(input_train_sol)), , drop = FALSE]

cat("Predicting mod1 DR of test cells\n")
pred_mod1 <- apply(dr_x1_train, 2, function(yi) {
  FNN::knn.reg(
    train = dr_x2_train,
    test = dr_x2_test,
    y = yi,
    k = min(15, nrow(dr_x1_test))
  )$pred
})

cat("Predicting mod2 DR of test cells\n")
pred_mod2 <- apply(dr_x2_train, 2, function(yi) {
  FNN::knn.reg(
    train = dr_x1_train,
    test = dr_x1_test,
    y = yi,
    k = min(15, nrow(dr_x1_test))
  )$pred
})

cat("Minimising distances between mod1 and mod2 pairs\n")
gen_vec <- function(z) {
  int <- seq_len(nrow(pred_mod1))

  i <- j <- c()
  resti <- int
  restj <- int

  while (length(resti) > 0) {
    ixi <- sample.int(length(resti), 1)
    newi <- resti[[ixi]]
    d1 <- proxy::dist(pred_mod1[restj, , drop = FALSE], dr_x1_test[newi, , drop = FALSE], method = "euclidean")
    d2 <- proxy::dist(pred_mod2[restj, , drop = FALSE], dr_x2_test[newi, , drop = FALSE], method = "euclidean")
    d12 <- d1 + d2
    ixj <- which.min(d12[, 1])
    newj <- restj[[ixj]]
    resti <- resti[-ixi]
    restj <- restj[-ixj]
    i <- c(i, newi)
    j <- c(j, newj)

    #  tibble(i, j); tibble(resti, restj)
  }

  tibble::tibble(i, j)
}

outs <- pbapply::pblapply(seq_len(par$n_pop), cl = n_cores, gen_vec)
# outs <- lapply(seq_len(par$n_pop), gen_vec)
df <- bind_rows(outs) %>%
  group_by(i, j) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n)) %>%
  mutate(gold = i == j)

knn_mat <- Matrix::sparseMatrix(
  i = df$i,
  j = df$j,
  x = df$n,
  dims = list(nrow(dr_x1_test), nrow(dr_x2_test))
)

# normalise to make rows sum to 1
rs <- Matrix::rowSums(knn_mat)
knn_mat@x <- knn_mat@x / rs[knn_mat@i + 1]

cat("Creating output anndata\n")
out <- anndata::AnnData(
  X = as(knn_mat, "CsparseMatrix"),
  uns = list(
    dataset_id = input_train_sol$uns[["dataset_id"]],
    method_id = method_id
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")
