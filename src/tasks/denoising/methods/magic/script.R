library(tidyverse)
library(reticulate)
library(anndata)

scprep <- import("scprep")
magic <- import("magic")
np <- import("numpy")

par <- list(
  input_train = "resources_test/denoising/pancreas/train.h5ad"
)

input_train <- anndata::read_h5ad(par$input_train)

X <- input_train$layers[["counts"]]

norm_out <- scprep$normalize$library_size_normalize(X, rescale = 10000, return_library_size = TRUE)

X_norm <- norm_out[[1]]
libsize <- norm_out[[2]] / 10000

# verify trans
qplot(
  as.vector(as.matrix(X)),
  as.vector(as.matrix(X_norm))
)

X_norm_unnorm <- scprep$utils$matrix_vector_elementwise_multiply(X_norm, libsize, axis = 0)
qplot(
  as.vector(as.matrix(X)),
  as.vector(as.matrix(X_norm_unnorm))
)

X_norm_trans <- scprep$utils$matrix_transform(X_norm, np$sqrt)

# verify
ix <- seq(min(X_norm), max(X_norm), by = .1)
patchwork::wrap_plots(
  qplot(
    as.vector(as.matrix(X_norm)),
    as.vector(as.matrix(X_norm_trans))
  ),
  qplot(ix, sqrt(ix)),
  ncol = 1
)

pheatmap::pheatmap(X_norm_trans)

# run magic
Y <- magic$MAGIC()$fit_transform(X_norm_trans, genes = "all_genes")

# verify
qplot(
  as.vector(as.matrix(X_norm_trans)),
  as.vector(as.matrix(Y))
)

# undo trans
Y_untrans <- scprep$utils$matrix_transform(Y, np$square)
qplot(
  as.vector(as.matrix(Y_untrans)),
  as.vector(as.matrix(Y))
)
qplot(
  as.vector(as.matrix(X_norm)),
  as.vector(as.matrix(Y_untrans))
)

Y_untrans_denorm <- scprep$utils$matrix_vector_elementwise_multiply(Y_untrans, libsize, axis = 0)
qplot(
  as.vector(as.matrix(X)),
  as.vector(as.matrix(Y_untrans_denorm))
)
