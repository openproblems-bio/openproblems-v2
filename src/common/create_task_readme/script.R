library(rlang)
library(purrr)
library(dplyr)

## VIASH START
par <- list(
  "task" = "batch_integration",
  "output" = "src/tasks/batch_integration/README.qmd",
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

# todo: fix hard coded node
order <- names(igraph::bfs(task_api$task_graph, "file_common_dataset")$order)

# render api details
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

# render authors
authors_str <- 
  if (nrow(task_api$authors) > 0) {
    paste0(
      "\n## Authors & contributors\n\n",
      task_api$authors %>% knitr::kable() %>% paste(collapse = "\n"),
      "\n"
    )
  } else {
    ""
  }

output <- strip_margin(glue::glue("
  §---
  §title: \"{task_api$task_info$label}\"
  §format: gfm
  §---
  §
  §{task_api$task_info$summary}
  §
  §Path: [`{task_dir}`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/{task_dir})
  §
  §## Motivation
  §
  §{task_api$task_info$motivation}
  §
  §## Description
  §
  §{task_api$task_info$description}
  §{authors_str}
  §## API
  §
  §{r_graph}
  §
  §{paste(r_details, collapse = '\n\n')}
  §
  §"), symbol = "§")

readr::write_lines(output, par$output)
