cat("Loading dependencies\n")
requireNamespace("anndata", quietly = TRUE)

## VIASH START
par <- list(
  input_test_sol = "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.train_sol.h5ad",
  output = "output.h5ad"
)
meta <- list(functionality_name = "foo")
## VIASH END

cat("Reading h5ad files\n")
input_test_sol <- anndata::read_h5ad(par$input_test_sol)

input_test_sol$uns[["method_id"]] <- meta$functionality_name

cat("Writing predictions to file\n")
zzz <- input_test_sol$write_h5ad(par$output, compression = "gzip")
