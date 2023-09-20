cat("Loading dependencies\n")
suppressPackageStartupMessages({
  requireNamespace("anndata", quietly = TRUE)
  library(Matrix, warn.conflicts = FALSE)
  requireNamespace("rliger", quietly = TRUE)
  library(SingleCellExperiment, warn.conflicts = FALSE)
  requireNamespace("Seurat", quietly = TRUE)
  requireNamespace("zellkonverter", quietly = FALSE)
})
## VIASH START
par <- list(
  input = 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
  output = 'output.h5ad',
  hvg = FALSE
)
meta <- list(
  functionality_name = "liger_embed"
)
## VIASH END

cat("Read input\n")
adata <- anndata::read_h5ad(par$input)
adata$X <- adata$layers[["normalized"]]
sce <- zellkonverter::AnnData2SCE(adata, skip_assays = FALSE, hdf5_backed = TRUE)

cat("Prepare data...\n")
sobj <- Seurat::as.Seurat(sce, data = NULL)
sobj@assays$RNA <- sobj@assays$originalexp
sobj@assays$RNA@counts <- sobj@assays$RNA@data


# Create Liger object
lobj <- rliger::seuratToLiger(
  sobj,
  combined.seurat = TRUE,
  meta.var = 'batch',
  raw.assay='RNA',
  renormalize = FALSE,
  remove.missing = FALSE
)

# We only pass nomarlized data, so store it as such
lobj@norm.data <- lobj@raw.data

cat("Run liger...\n")
# # Assign hvgs
# lobj@var.genes <- rownames(sobj@assays$RNA)

cat("scaleNotCenter\n")
lobj <- rliger::scaleNotCenter(
  lobj,
  remove.missing = FALSE,
)

cat("optimizeALS\n")
# Use tutorial coarse k suggests of 20.
print(lobj)
lobj <- rliger::optimizeALS(
  lobj,
  k = 20,
  thresh = 3,
  nrep = 5e-5,
  remove.missing = FALSE,
  verbose = TRUE
)

cat("quantileAlignSNF\n")
lobj <- rliger::quantileAlignSNF(
  lobj,
  resolution = 0.4,
  small.clust.thresh = 20,
  remove.missing = FALSE
)

cat("Reformat output\n")
adata$obsm[["X_emb"]] <- lobj@H.norm

cat("Store outputs...\n")
adata$uns[["method_id"]] <- meta$functionality_name
zzz <- adata$write_h5ad(par$output, compression = "gzip")
