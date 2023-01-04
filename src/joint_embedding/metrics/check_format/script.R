cat("Load dependencies\n")
library(assertthat, quietly = TRUE, warn.conflicts = FALSE)
library(anndata, warn.conflicts = FALSE)

## VIASH START
task <- "joint_embedding"
par <- list(
  input_solution = paste0("resources_test/", task, "/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.solution.h5ad"),
  input_prediction = paste0("resources_test/", task, "/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.prediction.h5ad"),
  output = paste0("resources_test/", task, "/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.scores.h5ad")
)
## VIASH END

cat("Read prediction h5ad\n")
ad_sol <- read_h5ad(par$input_solution)

cat("Checking solution h5ad\n")
correct_format <- tryCatch({
  # read prediction
  ad_pred <- read_h5ad(par$input_prediction)

  # check dataset id
  dataset_id <- ad_pred$uns[["dataset_id"]]
  assert_that(dataset_id == ad_sol$uns[["dataset_id"]])

  # check method id
  method_id <- ad_pred$uns[["method_id"]]
  assert_that(
    is.character(method_id),
    method_id != ""
  )

  # check X
  assert_that(
    ad_pred$n_obs == ad_sol$n_obs,
    ad_pred$n_vars >= 1,
    ad_pred$n_vars <= 100,
    !is.null(ad_pred$obs_names),
    all(ad_pred$obs_names == ad_sol$obs_names)
  )

  1
}, error = function(e) {
  cat("ERROR: ", e$message, "\n", sep = "")
  0
})


cat("Create output object\n")
out <- AnnData(
  shape = c(0, 0),
  uns = list(
    dataset_id = ad_pred$uns$dataset_id,
    method_id = ad_pred$uns$method_id,
    metric_ids = c("finished", "correct_format"),
    metric_values = c(1, correct_format)
  )
)

cat("Write output to h5ad file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")
