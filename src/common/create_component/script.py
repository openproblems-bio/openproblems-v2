from typing import Any
from pathlib import Path
import sys
import os
import re

## VIASH START
par = {
  "language": "python",
  "name": "new_comp",
  "output": "src/task/method/new_comp",
  "api_file": "src/task/api/comp_method.yaml",
  "viash_yaml": "_viash.yaml"
}
## VIASH END

# import helper function
sys.path.append(meta["resources_dir"])
from read_and_merge_yaml import read_and_merge_yaml

def strip_margin(text: str) -> str:
  return re.sub("(^|\n)[ \t]*\|", "\\1", text)

def create_config(par, component_type, pretty_name, script_path) -> str:
  info_str = generate_info(par, component_type, pretty_name)
  resources_str = generate_resources(par, script_path)
  docker_platform = generate_docker_platform(par)

  return strip_margin(f'''\
    |# The API specifies which type of component this is.
    |# It contains specifications for:
    |#   - The input/output files
    |#   - Common parameters
    |#   - A unit test
    |__merge__: {os.path.relpath(par["api_file"], par["output"])}
    |
    |functionality:
    |  # A unique identifier for your component (required).
    |  # Can contain only lowercase letters or underscores.
    |  name: {par["name"]}
    |
    |  # Metadata for your component
    |  info:
    |{info_str}
    |  # Component-specific parameters (optional)
    |  # arguments:
    |  #   - name: "--n_neighbors"
    |  #     type: "integer"
    |  #     default: 5
    |  #     description: Number of neighbors to use.
    |
    |  # Resources required to run the component
    |  resources:
    |{resources_str}
    |platforms:
    |  # Specifications for the Docker image for this component.
    |{docker_platform}
    |  # This platform allows running the component natively
    |  - type: native
    |  # Allows turning the component into a Nextflow module / pipeline.
    |  - type: nextflow
    |    directives:
    |      label: [midtime,midmem,midcpu]
    |'''
  )

def generate_info(par, component_type, pretty_name) -> str:
  """Generate the functionality info for a component."""
  if component_type in ["method", "control_method"]:
    str = strip_margin(f'''\
      |    # A relatively short label, used when rendering visualisations (required)
      |    label: {pretty_name}
      |    # A one sentence summary of how this method works (required). Used when 
      |    # rendering summary tables.
      |    summary: "FILL IN: A one sentence summary of this method."
      |    # A multi-line description of how this component works (required). Used
      |    # when rendering reference documentation.
      |    description: |
      |      FILL IN: A (multi-line) description of how this method works.
      |    # Which normalisation method this component prefers to use (required).
      |    preferred_normalization: log_cp10k
      |''')
    if component_type == "method":
      str += strip_margin(f'''\
        |    # A reference key from the bibtex library at src/common/library.bib (required).
        |    reference: bibtex_reference_key
        |    # URL to the documentation for this method (required).
        |    documentation_url: https://url.to/the/documentation
        |    # URL to the code repository for this method (required).
        |    repository_url: https://github.com/organisation/repository
        |''')
    return str
  elif component_type == "metric":
    return strip_margin(f'''\
      |    metrics:
      |      # A unique identifier for your metric (required).
      |      # Can contain only lowercase letters or underscores.
      |      name: {par["name"]}
      |      # A relatively short label, used when rendering visualisarions (required)
      |      label: {pretty_name}
      |      # A one sentence summary of how this metric works (required). Used when 
      |      # rendering summary tables.
      |      summary: "FILL IN: A one sentence summary of this metric."
      |      # A multi-line description of how this component works (required). Used
      |      # when rendering reference documentation.
      |      description: |
      |        FILL IN: A (multi-line) description of how this metric works.
      |      # A reference key from the bibtex library at src/common/library.bib (required).
      |      reference: bibtex_reference_key
      |      # URL to the documentation for this metric (required).
      |      documentation_url: https://url.to/the/documentation
      |      # URL to the code repository for this metric (required).
      |      repository_url: https://github.com/organisation/repository
      |      # The minimum possible value for this metric (required)
      |      min: 0
      |      # The maximum possible value for this metric (required)
      |      max: 1
      |      # Whether a higher value represents a 'better' solution (required)
      |      maximize: true
      |''')


def generate_resources(par, script_path) -> str:
  """Add the script to the functionality resources."""
  if par["language"] == "python":
    type_str = "python_script"
  elif par["language"] == "r":
    type_str = "r_script"

  return strip_margin(f'''\
    |    # The script of your component (required)
    |    - type: {type_str}
    |      path: {script_path}
    |    # Additional resources your script needs (optional)
    |    # - type: file
    |    #   path: weights.pt
    |''')

def generate_docker_platform(par) -> str:
  """Set up the docker platform for Python."""
  if par["language"] == "python":
    image_str = "ghcr.io/openproblems-bio/base_python:1.0.4"
    extra = ""
  elif par["language"] == "r":
    image_str = "ghcr.io/openproblems-bio/base_r:1.0.4"
    extra = strip_margin(f'''\
      |      - type: r
      |        packages: [ arrow, readr ]
      |''')
  return strip_margin(f'''\
    |  - type: docker
    |    image: {image_str}
    |    # Add custom dependencies here (optional). For more information, see
    |    # https://viash.io/reference/config/platforms/docker/#setup .
    |    setup:
    |      - type: python
    |        packages: [ fastparquet ]
    |{extra}''')

def set_par_values(config) -> None:
  """Adds values to each of the arguments in a config file."""
  args = config['functionality']['arguments']
  for argi, arg in enumerate(args):
    key = re.sub("^-*", "", arg['name'])

    # find value
    if arg["type"] != "file":
      value = arg.get("default", arg.get("example", "..."))
    elif key == "de_train":
      value = "resources/neurips-2023-kaggle/de_train.parquet"
    elif key == "de_train_h5ad":
      value = "resources/neurips-2023-kaggle/2023-09-12_de_by_cell_type_train.h5ad"
    elif key == "id_map":
      value = "resources/neurips-2023-kaggle/id_map.csv"
    else:
      key_strip = key.replace("output_", "")
      value = f'{key_strip}.h5ad'

    # store key and value
    config['functionality']['arguments'][argi]["key"] = key
    config['functionality']['arguments'][argi]["value"] = value
  

def create_python_script(par, config, type):
  script = strip_margin('''\
    |import pandas as pd
    |
    |## VIASH START
    |par = {
    |  "de_train": "resources/neurips-2023-kaggle/de_train.parquet",
    |  "de_test": "resources/neurips-2023-kaggle/de_test.parquet",
    |  "id_map": "resources/neurips-2023-kaggle/id_map.csv",
    |  "output": "output.parquet",
    |}
    |## VIASH END
    |
    |print('Reading input files', flush=True)
    |de_train = pd.read_parquet(par["de_train"])
    |id_map = pd.read_csv(par["id_map"])
    |gene_names = [col for col in de_train.columns if col not in {"cell_type", "sm_name", "sm_lincs_id", "SMILES", "split", "control", "index"}]
    |
    |print('Preprocess data', flush=True)
    |# ... preprocessing ...
    |
    |print('Train model', flush=True)
    |# ... train model ...
    |
    |print('Generate predictions', flush=True)
    |# ... generate predictions ...
    |
    |print('Write output to file', flush=True)
    |output = pd.DataFrame(
    |  # ... TODO: fill in data ...
    |  index=id_map["id"],
    |  columns=gene_names
    |).reset_index()
    |output.to_parquet(par["output"])
    |''')

  return script

def create_r_script(par, api_spec, type):
  script = strip_margin(f'''\
    |requireNamespace("arrow", quietly = TRUE)
    |requireNamespace("readr", quietly = TRUE)
    |
    |## VIASH START
    |par <- list(
    |  de_train = "resources/neurips-2023-kaggle/de_train.parquet",
    |  id_map = "resources/neurips-2023-kaggle/id_map.csv",
    |  output = "output.parquet"
    |)
    |## VIASH END
    |
    |cat("Reading input files\\n")
    |de_train <- arrow::read_parquet(par$de_train)
    |id_map <- readr::read_csv(par$id_map)
    |
    |cat("Preprocess data\\n")
    |# ... preprocessing ...
    |
    |cat("Train model\\n")
    |# ... train model ...
    |
    |cat("Generate predictions\\n")
    |# ... generate predictions ...
    |
    |cat("Write output to file\\n")
    |output <- data.frame(
    |  id = id_map$id,
    |  # ... more columns ...
    |  check.names = FALSE
    |)
    |arrow::write_parquet(output, par$output)
    |''')

  return script



def main(par):
  ####### CHECK INPUTS #######
  print("Check inputs", flush=True)
  assert re.match("[a-z][a-z0-9_]*", par["name"]), "Name should match the regular expression '[a-z][a-z0-9_]*'. Example: 'my_component'."
  assert len(par['name']) <= 50, "Method name should be at most 50 characters."

  pretty_name = re.sub("_", " ", par['name']).title()

  ####### CHECK LANGUAGE #######
  print("Check language", flush=True)
  # check language and determine script path
  if par["language"] == "python":
    script_path = "script.py"
  elif par["language"] == "r":
    script_path = "script.R"
  else:
    sys.exit(f"Unrecognized language parameter '{par['language']}'.")

  ## CHECK API FILE
  print("Check API file", flush=True)
  api_file = Path(par["api_file"])
  viash_yaml = Path(par["viash_yaml"])
  project_dir = viash_yaml.parent
  if not api_file.exists():
    comp_types = [x.with_suffix("").name.removeprefix("comp_") for x in api_file.parent.glob("**/comp_*.y*ml")]
    list.sort(comp_types)
    sys.exit(strip_margin(f"""\
      |Error: Invalid --type argument.
      |  Reason: Could not find API file at '{api_file.relative_to(project_dir)}'.
      |  Possible values for --type: {', '.join(comp_types)}."""))
  
  ## READ API FILE
  print("Read API file", flush=True)
  api = read_and_merge_yaml(api_file)
  comp_type = api.get("functionality", {}).get("info", {}).get("type", {})
  if not comp_type:
    sys.exit(strip_margin(f"""\
      |Error: API file is incorrectly formatted.
      |  Reason: Could not find component type at `.functionality.info.type`.'
      |  Please fix the formatting of the API file."""))

  ####### CREATE OUTPUT DIR #######
  print("Create output dir", flush=True)
  out_dir = Path(par["output"])
  out_dir.mkdir(exist_ok=True)

  ####### CREATE CONFIG #######
  print("Create config", flush=True)
  config_file = out_dir / "config.vsh.yaml"

  # get config template
  config_str = create_config(par, comp_type, pretty_name, script_path)

  with open(config_file, "w") as f:
    f.write(config_str)

  ####### CREATE SCRIPT #######
  print("Create script", flush=True)
  script_file = out_dir / script_path

  # set reasonable values
  set_par_values(api)

  if par["language"] == "python":
    script_out = create_python_script(par, api, comp_type)

  if par["language"] == "r":
    script_out = create_r_script(par, api, comp_type)
  
  # write script
  with open(script_file, "w") as f:
    f.write(script_out)

  print("Done!", flush=True)


if __name__ == "__main__":
  main(par)