library(testthat, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)

par <- list(
  input_mod1 = "resources_test/common/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.output_rna.h5ad",
  input_mod2 = "resources_test/common/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.output_mod2.h5ad",
  output_mod1 = "output_mod1.h5ad",
  output_mod2 = "output_mod2.h5ad",
  output_solution = "solution.h5ad"
)

cat("> Running censor component\n")
out <- processx::run(
  command = paste0("./", meta["functionality_name"]),
  args = c(
    "--input_mod1", par$input_mod1,
    "--input_mod2", par$input_mod2,
    "--output_mod1", par$output_mod1,
    "--output_mod2", par$output_mod2,
    "--output_solution", par$output_solution
  ),
  stderr_to_stdout = TRUE
)

cat("> Checking whether output files were created\n")
expect_true(file.exists(par$output_mod1))
expect_true(file.exists(par$output_mod2))
expect_true(file.exists(par$output_solution))

cat("> Reading h5ad files\n")
input_mod1 <- anndata::read_h5ad(par$input_mod1)
input_mod2 <- anndata::read_h5ad(par$input_mod2)
output_mod1 <- anndata::read_h5ad(par$output_mod1)
output_mod2 <- anndata::read_h5ad(par$output_mod2)
output_solution <- anndata::read_h5ad(par$output_solution)

cat("> Checking contents of h5ad files\n")
expect_equal(output_mod1$uns[["dataset_id"]], paste0(input_mod1$uns[["dataset_id"]], "_JE"))
expect_equal(output_mod2$uns[["dataset_id"]], paste0(input_mod1$uns[["dataset_id"]], "_JE"))
expect_equal(output_solution$uns[["dataset_id"]], paste0(input_mod1$uns[["dataset_id"]], "_JE"))
expect_equal(output_mod1$uns[["organism"]], input_mod1$uns[["organism"]])
expect_equal(output_mod2$uns[["organism"]], input_mod1$uns[["organism"]])
expect_equal(output_solution$uns[["organism"]], input_mod1$uns[["organism"]])
expect_equal(output_mod1$n_obs, input_mod1$n_obs)
expect_equal(output_mod2$n_obs, input_mod2$n_obs)
expect_equal(output_mod1$n_vars, input_mod1$n_vars)
expect_equal(output_mod2$n_vars, input_mod2$n_vars)
expect_equal(output_mod1$var_names, input_mod1$var_names)
expect_equal(output_mod2$var_names, input_mod2$var_names)
expect_equal(output_mod1$obs_names, input_mod1$obs_names)
expect_equal(output_mod2$obs_names, input_mod2$obs_names)

# TODO check contents of matrices, check rownames

cat("> Test succeeded!\n")
