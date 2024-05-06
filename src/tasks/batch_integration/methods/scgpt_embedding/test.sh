#!/bin/bash

set -e

viash run src/tasks/batch_integration/methods/scgpt/config.vsh.yaml -- \
  --input "resources_test/batch_integration/pancreas/dataset.h5ad" \
  --model "resources_test/scGPT_human/best_model.pt" \
  --model_config "resources_test/scGPT_human/args.json" \
  --model_vocab "resources_test/scGPT_human/vocab.json" \
  --output "output/temp/scgpt/pancreas.h5ad" \