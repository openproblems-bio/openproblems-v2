library(RcppDist)

source("opt/BOOST-GP/boost.gp.R")

# VIASH START
par <- list(
    "input_data" = "resources_test/spatially_variable_genes/10x_visium_mouse_brain/input_data.h5ad",
    "output" = "output.h5ad"
)
meta <- list(
    "functionality_name" = "BOOST-GP"
)
# VIASH END

cat("Load data")
counts <- as.matrix(adata$layers[["counts"]])
colnames(counts) <- adata$var_names
rownames(counts) <- adata$obs_names
mode(counts) <- "integer"

loc <- as.data.frame(adata$obsm[["spatial"]])
rownames(loc) <- adata$obs_names
colnames(loc) <- c("x", "y")

cat("Run BOOST-GP")
df <- as.data.frame(boost.gp(Y = counts, loc = loc, iter = 10, burn = 5))

df$feature_name <- rownames(df)

df <- subset(df, select = c("feature_name", "PPI"))
colnames(df) <- c("feature_name", "pred_spatial_var_score")

# save output
cat("Write output AnnData to file\n")
output <- anndata::AnnData(
    shape = adata$shape,
    var = df,
    uns = list(
        "dataset_id" = adata$uns["dataset_id"],
        "method_id" = meta["functionality_name"]
    )
)

anndata::write_h5ad(anndata = output, filename = par$output)
