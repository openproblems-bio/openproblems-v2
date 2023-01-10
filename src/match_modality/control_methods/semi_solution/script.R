cat("Loading dependencies\n")
requireNamespace("anndata", quietly = TRUE)

## VIASH START
par <- list(
  input_test_sol = "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_sol.h5ad",
  output = "output.h5ad"
)
meta <- list(functionality_name = "foo")
## VIASH END

cat("Reading h5ad files\n")
input_test_sol <- anndata::read_h5ad(par$input_test_sol)

# randomly fill in gold standard values
input_test_sol$X@x <- runif(length(input_test_sol$X@x))

# fill other values with random values as well
ix <- sample.int(nrow(input_test_sol) * ncol(input_test_sol), nrow(input_test_sol) * 10)
input_test_sol$X[ix] <- runif(length(ix))

# normalise to make rows sum to 1
rs <- Matrix::rowSums(input_test_sol$X)
input_test_sol$X@x <- input_test_sol$X@x / rs[input_test_sol$X@i + 1]

input_test_sol$uns[["method_id"]] <- meta$functionality_name

cat("Writing predictions to file\n")
zzz <- input_test_sol$write_h5ad(par$output, compression = "gzip")
