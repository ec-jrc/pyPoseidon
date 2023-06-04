#!/usr/bin/env bash
#

set -xeuo pipefail

# Binary
./scripts/merge_environment_yaml.py ./dependencies/python3.9.yml   ./dependencies/binary.yml > ./environments/binary-p3.9.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.10.yml  ./dependencies/binary.yml > ./environments/binary-p3.10.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.11.yml  ./dependencies/binary.yml > ./environments/binary-p3.11.yml

# Base
./scripts/merge_environment_yaml.py ./dependencies/python3.9.yml  ./dependencies/binary.yml ./dependencies/main.yml > ./environments/base-p3.9.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.10.yml ./dependencies/binary.yml ./dependencies/main.yml > ./environments/base-p3.10.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.11.yml ./dependencies/binary.yml ./dependencies/main.yml > ./environments/base-p3.11.yml

# Viz
./scripts/merge_environment_yaml.py ./dependencies/python3.9.yml  ./dependencies/binary.yml ./dependencies/main.yml ./dependencies/viz.yml > ./environments/viz-p3.9.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.10.yml ./dependencies/binary.yml ./dependencies/main.yml ./dependencies/viz.yml > ./environments/viz-p3.10.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.11.yml ./dependencies/binary.yml ./dependencies/main.yml ./dependencies/viz.yml > ./environments/viz-p3.11.yml

# Full
./scripts/merge_environment_yaml.py ./dependencies/python3.9.yml  ./dependencies/binary.yml ./dependencies/dev.yml ./dependencies/docs.yml ./dependencies/main.yml ./dependencies/viz.yml > ./environments/full-p3.9.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.10.yml ./dependencies/binary.yml ./dependencies/dev.yml ./dependencies/docs.yml ./dependencies/main.yml ./dependencies/viz.yml > ./environments/full-p3.10.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.11.yml ./dependencies/binary.yml ./dependencies/dev.yml ./dependencies/docs.yml ./dependencies/main.yml ./dependencies/viz.yml > ./environments/full-p3.11.yml
