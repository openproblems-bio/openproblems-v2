library(rlang)
library(purrr)
library(dplyr)

## VIASH START
par <- list(
  "task" = "denoising",
  "output" = "src/tasks/denoising/README.qmd",
  "viash_yaml" = "_viash.yaml"
)
meta <- list(
  "resources_dir" = "src/common/helper_functions"
)
## VIASH END

# import helper function
source(paste0(meta["resources_dir"], "/read_and_merge_yaml.R"))
source(paste0(meta["resources_dir"], "/strip_margin.R"))
source(paste0(meta["resources_dir"], "/read_api_files.R"))

# find task dir
task_dir <- paste0(dirname(par[["viash_yaml"]]), "/src/tasks/", par[["task"]]) %>%
  gsub("^\\./", "", .)
task_api <- read_task_api(task_dir)


r_graph <- render_task_graph(task_api)

order <- names(igraph::bfs(task_api$task_graph, "anndata_common_dataset")$order)

r_details <- map_chr(
  order,
  function(file_name) {
    if (file_name %in% names(task_api$comp_specs)) {
      render_component(task_api$comp_specs[[file_name]])
    } else {
      render_file(task_api$file_specs[[file_name]])
    }
  }
)

## TODO: render authors

output <- strip_margin(glue::glue("
  §---
  §title: \"{task_api$task_info$task_name}\"
  §format: gfm
  §---
  §
  §Path: [`{task_dir}`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/{task_dir})
  §
  §{task_api$task_info$description}
  §
  §{r_graph}
  §
  §{paste(r_details, collapse = '\n\n')}
  §
  §"), symbol = "§")

readr::write_lines(output, par$output)
