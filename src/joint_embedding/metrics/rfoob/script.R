cat("Load dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(testthat, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)

## VIASH START
par <- list(
  input_solution = "resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.solution.h5ad",
  input_prediction = "resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.prediction.h5ad",
  output = "scores.h5ad"
)
## VIASH END

cat("Read solution h5ad\n")
ad_sol <- anndata::read_h5ad(par$input_solution)

cat("Read prediction h5ad\n")
expect_true(
  grepl("\\.h5ad$", par$input_prediction),
  info = "Prediction file should be an h5ad file"
)
ad_pred <-
  tryCatch({
    anndata::read_h5ad(par$input_prediction)
  }, error = function(e) {
    stop(paste0("Can't open prediction h5ad file. Detailed error message:\n", e$message))
  })
expect_true(
  ad_sol$uns$dataset_id == ad_pred$uns$dataset_id
)

cat("Calculating metrics\n")
df <- data.frame(as.matrix(ad_pred$obsm[["X_emb"]]), SOLUTION_CELL_TYPE = ad_sol$obs[["cell_type"]])
rf1 <- ranger::ranger(SOLUTION_CELL_TYPE ~ ., df)

df <- data.frame(as.matrix(ad_pred$obsm[["X_emb"]]), SOLUTION_PSEUDOTIME_ORDER = ad_sol$obs$pseudotime_order_GEX)
df <- df[is.finite(df$SOLUTION_PSEUDOTIME_ORDER), , drop = FALSE]
rf2 <- ranger::ranger(SOLUTION_PSEUDOTIME_ORDER ~ ., df)

colname <- colnames(ad_sol$obs)[grepl("pseudotime_order_A.*", colnames(ad_sol$obs))]
df <- data.frame(as.matrix(ad_pred$obsm[["X_emb"]]), SOLUTION_PSEUDOTIME_ORDER = ad_sol$obs[[colname]])
df <- df[is.finite(df$SOLUTION_PSEUDOTIME_ORDER), , drop = FALSE]
rf3 <- ranger::ranger(SOLUTION_PSEUDOTIME_ORDER ~ ., df)

df <- data.frame(as.matrix(ad_pred$obsm[["X_emb"]]), SOLUTION_BATCH = ad_sol$obs$batch)
rf4 <- ranger::ranger(SOLUTION_BATCH ~ ., df)

metric_values <- c(
  rfoob_celltype_accuracy = 1 - rf1$prediction.error,
  rfoob_pseudotimegex_rsq = rf2$r.squared,
  rfoob_pseudotimeadt_rsq = rf3$r.squared,
  rfoob_batch_error = rf4$prediction.error
)

cat("Create output object\n")
out <- anndata::AnnData(
  shape = c(0, 0),
  uns = list(
    dataset_id = ad_pred$uns$dataset_id,
    method_id = ad_pred$uns$method_id,
    metric_ids = names(metric_values),
    metric_values = metric_values
  )
)

cat("Write output to h5ad file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")
