library(assertthat, quietly = TRUE)
library(Matrix, quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)

## VIASH START
# This code block will be replaced by viash at runtime.
meta <- list(functionality_name = "foo")
## VIASH END

# determine filenames and arguments
testpar <- list(
  "input_solution" = "temp_sol.h5ad",
  "input_prediction" = "temp_pred.h5ad",
  "output" = "temp_out.h5ad"
)
command <- paste0("./", meta[["functionality_name"]])
args <- unlist(rbind(paste0("--", names(testpar)), unname(testpar)))

# uncomment this for manual testing
# command <- "viash"
# args <- c("run", "src/match_modality/metrics/aupr/config.vsh.yaml", "--", args)

cat("Creating test files\n")
ad_sol <- anndata::AnnData(
  X = as(Matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, byrow = TRUE, sparse = TRUE), "CsparseMatrix"),
  layers = list(
    neighbors = as(Matrix(c(1, 0.5, 0, 0.5, 1, 0.25, 0, 0.25, 1), nrow = 3, byrow = TRUE, sparse = TRUE), "CsparseMatrix")
  ),
  uns = list(dataset_id = "simple", pairing_ix = c(0, 1, 2)),
  obs = data.frame(
    cell_type = c("a", "a", "b")
  )
)
ad_pred <- anndata::AnnData(
  X = as(Matrix(c(1, .1, .2, .3, .9, .4, .5, .6, .8), nrow = 3, byrow = TRUE, sparse = TRUE), "CsparseMatrix"),
  uns = list(dataset_id = "simple", method_id = "simple")
)

ad_sol$write_h5ad(testpar$input_solution, compression = "gzip")
ad_pred$write_h5ad(testpar$input_prediction, compression = "gzip")

cat("> Running metric\n")
out <- processx::run(
  command = command,
  args = args,
  stderr_to_stdout = TRUE
)

cat("> Reading metric scores\n")
assert_that(file.exists(testpar$output))
ad_out <- anndata::read_h5ad(testpar$output)

scores1 <- ad_out$uns$metric_values
names(scores1) <- ad_out$uns$metric_ids
# assert_that(
#   scores1[["pairing_aupr"]] >= scores1[["neighbor_aupr"]],
#   scores1[["neighbor_aupr"]] >= scores1[["celltype_aupr"]],
#   scores1[["pairing_auroc"]] >= scores1[["neighbor_auroc"]],
#   scores1[["pairing_auroc"]] >= scores1[["celltype_auroc"]]
# )


cat("Creating test files\n")
# pairing_ix <- c(2, 1, 3)
pairing_ix <- c(3, 1, 2)
ad_sol <- anndata::AnnData(
  X = as(diag(length(pairing_ix)), "CsparseMatrix")[,pairing_ix],
  layers = list(
    neighbors = as(Matrix(c(1, 0.5, 0, 0.5, 1, 0.25, 0, 0.25, 1), nrow = 3, byrow = TRUE, sparse = TRUE), "CsparseMatrix")[,pairing_ix]
  ),
  uns = list(dataset_id = "simple", pairing_ix = pairing_ix-1),
  obs = data.frame(
    cell_type = c("a", "a", "b")
  )
)
ad_pred <- anndata::AnnData(
  X = as(Matrix(c(1, .1, .2, .3, .9, .4, .5, .6, .8), nrow = 3, byrow = TRUE, sparse = TRUE), "CsparseMatrix")[,pairing_ix],
  uns = list(dataset_id = "simple", method_id = "simple")
)

ad_sol$write_h5ad(testpar$input_solution, compression = "gzip")
ad_pred$write_h5ad(testpar$input_prediction, compression = "gzip")

cat("> Running metric\n")
out <- processx::run(
  command = command,
  args = args,
  stderr_to_stdout = TRUE
)

cat("> Reading metric scores\n")
assert_that(file.exists(testpar$output))
ad_out <- anndata::read_h5ad(testpar$output)

scores2 <- ad_out$uns$metric_values
names(scores2) <- ad_out$uns$metric_ids

assert_that(all(scores1 == scores2))

cat("> Test succeeded!\n")
