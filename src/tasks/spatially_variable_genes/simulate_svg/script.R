# suppressMessages(library(scDesign3))
# suppressMessages(library(scales))
# suppressMessages(library(dplyr))
# suppressMessages(library(ggplot2))
# suppressMessages(library(cowplot))
# suppressMessages(library(Seurat))
# suppressMessages(library(SingleCellExperiment))
requireNamespace("scDesign3", quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)
requireNamespace("Matrix", quietly = TRUE)
requireNamespace("SingleCellExperiment", quietly = TRUE)
library(rlang)

## VIASH START
par <- list(
  input = "dataset_raw.h5ad",
  output = "dataset_sim.h5ad"
)
meta <- list(
  cpus = 30L
)
## VIASH END

adata <- anndata::read_h5ad(par$input)

adata <- adata[,1:100]

# get data
df_loc <- as.data.frame(adata$obsm[['spatial']])
colnames(df_loc) <- c("spatial1", "spatial2")
rownames(df_loc) <- adata$obs_names

# transform into SCE
ref_sce <- SingleCellExperiment(
  list(counts = Matrix::t(adata$X)),
  colData = df_loc
)

ref_sce

# constructs the input data for fit_marginal.
ref_data <- scDesign3::construct_data(
  sce = ref_sce,
  assay_use = "counts",
  celltype = NULL,
  pseudotime = NULL,
  spatial = c("spatial1", "spatial2"),
  other_covariates = NULL,
  corr_by = "1"
)

# fit expression of each gene with GP model
ref_marginal <- scDesign3::fit_marginal(
  data = ref_data,
  predictor = "gene",
  mu_formula = "s(spatial1, spatial2, bs = 'gp', k = 500)",
  sigma_formula = "1",
  family_use = "nb",
  parallelization = "pbmcmapply",
  n_cores = meta$cpus %||% 2L,
  usebam = FALSE,
  trace = TRUE
)

# 
ref_copula <- scDesign3::fit_copula(
  sce = ref_sce,
  assay_use = "counts",
  marginal_list = ref_marginal,
  family_use = "nb",
  copula = "gaussian",
  parallelization = "pbmcmapply",
  n_cores = meta$cpus %||% 2L,
  input_data = ref_data$dat
)

ref_para <- scDesign3::extract_para(
  sce = ref_sce,
  marginal_list = ref_marginal,
  family_use = "nb",
  new_covariate = ref_data$newCovariate,
  data = ref_data$dat,
  parallelization = "pbmcmapply",
  n_cores = meta$cpus %||% 2L
)

dev_explain <- sapply(ref_marginal, function(x) {
  sum = summary(x$fit)
  return(sum$dev.expl)
})
dev_ordered <- order(dev_explain, decreasing = TRUE)
num_de <- 50
ordered <- dev_explain[dev_ordered]
sel_genes <- names(ordered)[1:num_de]

ref_sce <- ref_sce[sel_genes, ]

# constructs the input data for fit_marginal.
ref_data <- construct_data(
  sce = ref_sce,
  assay_use = "counts",
  celltype = NULL,
  pseudotime = NULL,
  spatial = c("spatial1", "spatial2"),
  other_covariates = NULL,
  corr_by = "1"
)

# fit expression of each gene with GP model
ref_marginal <- fit_marginal(
  data = ref_data,
  predictor = "gene",
  mu_formula = "s(spatial1, spatial2, bs = 'gp', k = 500)", 
  sigma_formula = "1",
  family_use = "nb",
  n_cores = 2,
  usebam = FALSE,
  trace = TRUE
)

ref_copula <- fit_copula(
  sce = ref_sce,
  assay_use = "counts",
  marginal_list = ref_marginal,
  family_use = "nb",
  copula = "gaussian",
  n_cores = 2,
  input_data = ref_data$dat
)

ref_para <- extract_para(
  sce = ref_sce,
  marginal_list = ref_marginal,
  n_cores = 5,
  family_use = "nb",
  new_covariate = ref_data$newCovariate,
  data = ref_data$dat
)

sim_count <- simu_new(
  sce = ref_sce,
  mean_mat = ref_para$mean_mat,
  sigma_mat = ref_para$sigma_mat,
  zero_mat = ref_para$zero_mat,
  quantile_mat = NULL,
  copula_list = ref_copula$copula_list,
  n_cores = 5,
  family_use = "nb",
  input_data = ref_data$dat,
  new_covariate = ref_data$newCovariate,
  important_feature = rep(TRUE, dim(ref_sce)[1]),
  filtered_gene = NULL
)

sim_sce <- SingleCellExperiment(list(counts = sim_count),
                                colData = ref_data$newCovariate)

# generate non-spatially variable mean values with shuffling
shuffle_idx <- sample(nrow(ref_para$mean_mat))
non_de_mat <- ref_para$mean_mat[shuffle_idx, ]

# simulate data with varied spatial variability
count <- lapply(seq(0, 1.0, 0.05), function(alpha){
    sim_count <- simu_new(
        sce = ref_sce,
        mean_mat = alpha * ref_para$mean_mat + (1 - alpha) * non_de_mat,
        sigma_mat = ref_para$sigma_mat,
        zero_mat = ref_para$zero_mat,
        quantile_mat = NULL,
        copula_list = ref_copula$copula_list,
        n_cores = 5,
        family_use = "nb",
        input_data = ref_data$dat,
        new_covariate = ref_data$newCovariate,
        important_feature = rep(TRUE, dim(ref_sce)[1]),
        filtered_gene = NULL)

    rownames(sim_count) <- paste0(rownames(sim_count), "_", alpha)
    
    return(sim_count)

}) %>% do.call(rbind, .)