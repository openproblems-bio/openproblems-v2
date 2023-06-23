#!/bin/bash

viash config view src/tasks/label_projection/workflows/run/config.vsh.yaml \
  --parse_argument_groups > \
  src/tasks/label_projection/workflows/run/.config.vsh.yaml