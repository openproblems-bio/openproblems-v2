requireNamespace("anndata", quietly = TRUE)
requireNamespace("SIMLR", quietly = TRUE)

## VIASH START
par <- list(
  input = "resources_test/dimensionality_reduction/pancreas/dataset.h5ad",
  output = "output.h5ad",
  n_clusters = NULL,
  n_dim = NULL,
  tuning_param = 10,
  impute = FALSE,
  normalize = FALSE,
  cores_ratio = 1
)
meta <- list(
  functionality_name = "simlr"
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

cat("Clean up memory\n")
normalized <- as.matrix(input$layers[["normalized"]])

output <- anndata::AnnData(
  uns = list(
    dataset_id = input$uns[["dataset_id"]],
    method_id = meta$functionality_name,
    normalization_id = input$uns[["normalization_id"]]
  ),
  shape = input$shape
)
rm(input)
gc()

if (is.null(par$n_clusters)) {
  cat("Estimating the number of clusters\n")
  set.seed(1)
  NUMC <- 2:5
  estimates <- SIMLR::SIMLR_Estimate_Number_of_Clusters(
    X = normalized,
    NUMC = NUMC,
    cores.ratio = par$cores_ratio
  )
  n_clusters <- NUMC[which.min(estimates$K2)]
} else {
  n_clusters <- par$n_clusters
}

if (is.null(par$n_dim)) {
  n_dim <- NA
} else {
  n_dim <- par$n_dim
}

cat("Running SIMLR\n")
simlr_result <- SIMLR::SIMLR(
  X = normalized,
  c = n_clusters,
  no.dim = n_dim,
  k = par$tuning_param,
  if.impute = par$impute,
  normalize = par$normalize,
  cores.ratio = par$cores_ratio
)

cat("Write output AnnData to file\n")
output$obsm$X_emb <- simlr_result$ydata
output$write_h5ad(par$output, compression = "gzip")
