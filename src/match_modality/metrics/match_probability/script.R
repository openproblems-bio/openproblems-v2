cat("Load dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(testthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)

## VIASH START
par <- list(
  input_solution = "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_sol.h5ad",
  input_prediction = "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.prediction.h5ad",
  output = "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.scores.h5ad"
)
## VIASH END

cat("Read solution h5ad\n")
ad_sol <- anndata::read_h5ad(par$input_solution, backed = "r")

cat("Read prediction h5ad\n")
ad_pred <- anndata::read_h5ad(par$input_prediction)

cat("Unscrambling predictions\n")
pairing_ix <- ad_sol$uns[["pairing_ix"]]
X_pred <- as(ad_pred$X, "CsparseMatrix")[, order(pairing_ix)]
dimnames(X_pred) <- list(NULL, NULL)

# set negative values to 0
X_pred@x <- ifelse(X_pred@x < 0, 0, X_pred@x)

cat("Calculating normalisation factors\n")
rowSum <- Matrix::rowSums(X_pred)

cat("Computing the match modality score\n")
match_probability_vec <- diag(X_pred) / rowSum

match_probability <- mean(match_probability_vec)

cat("Create output object\n")
out <- anndata::AnnData(
  shape = c(0, 0),
  uns = list(
    dataset_id = ad_pred$uns$dataset_id,
    method_id = ad_pred$uns$method_id,
    metric_ids = list("match_probability"),
    metric_values = list(match_probability),
    per_cell = list(
      match_probability = match_probability_vec
    )
  )
)

# should we also save the metrics object?
# this would allow for plotting the auroc and aupr curves afterwards.

cat("Write output to h5ad file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")