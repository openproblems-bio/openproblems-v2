library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(tibble, warn.conflicts = FALSE)
library(ontologyIndex, warn.conflicts = FALSE)

# february 2023 release is used
# https://github.com/obophenotype/cell-ontology/releases/download/v2023-02-15/cl.obo

## VIASH START
par <- list(
  input = "dataset_cxg2.h5ad",
  ontology = "cl.obo",
  input_term = "cell_type_ontology_term_id",
  struct = "obs",
  output = "output.h5ad",
  output_term = "cell_type_ontology_term_id",
  output_name = "cell_type",
  output_obsolete = "cell_type_ontology_obsolete"
)
## VIASH END

cat("Read ontology\n")
ont <- ontologyIndex::get_ontology(
  par$ontology,
  extract_tags = "everything"
)
ont_tib <- ont %>%
  as.data.frame %>%
  select(id, name, obsolete, replaced_by) %>%
  as_tibble

cat("Read anndata\n")
adata <- anndata::read_h5ad(par$input, backed = "r")

cat("Find terms\n")
term_ids <- adata[[par$struct]][[par$input_term]]

unique_term_ids <- as.character(unique(term_ids))

cat("Look for obsolete or replaced terms\n")
ont_map <- ont_tib %>%
  slice(match(unique_term_ids, id)) %>%
  transmute(orig_id = id, id = ifelse(replaced_by != "", replaced_by, id)) %>%
  left_join(ont_tib %>% select(id, name, obsolete), by = "id")

cat("Store new columns in data structure\n")
new_data <- ont_map %>% slice(match(term_ids, orig_id))
adata[[par$struct]][[par$output_term]] <- new_data$id
adata[[par$struct]][[par$output_name]] <- new_data$name
adata[[par$struct]][[par$output_obsolete]] <- new_data$obsolete

cat("Write to file\n")
# anndata::write_h5ad(adata, par$output, compression = "gzip")
anndata::write_h5ad(adata, par$output)
