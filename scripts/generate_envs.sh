#!/usr/bin/env bash
#

set -xeuo pipefail

# Binary
./scripts/merge_environment_yaml.py ./dependencies/python3.9.yml   ./dependencies/binary.yml | sed '/^$/d' > ./environments/binary-p3.9.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.10.yml  ./dependencies/binary.yml | sed '/^$/d' > ./environments/binary-p3.10.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.11.yml  ./dependencies/binary.yml | sed '/^$/d' > ./environments/binary-p3.11.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.12.yml  ./dependencies/binary.yml | sed '/^$/d' > ./environments/binary-p3.12.yml

# Binary + Telemac
./scripts/merge_environment_yaml.py ./dependencies/python3.9.yml   ./dependencies/binary.yml ./dependencies/telemac.yml > ./environments/binary-telemac-p3.9.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.10.yml  ./dependencies/binary.yml ./dependencies/telemac.yml > ./environments/binary-telemac-p3.10.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.11.yml  ./dependencies/binary.yml ./dependencies/telemac.yml > ./environments/binary-telemac-p3.11.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.12.yml  ./dependencies/binary.yml ./dependencies/telemac.yml > ./environments/binary-telemac-p3.12.yml

# Base
./scripts/merge_environment_yaml.py ./dependencies/python3.9.yml  ./dependencies/binary.yml ./dependencies/main.yml | sed '/^$/d' > ./environments/base-p3.9.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.10.yml ./dependencies/binary.yml ./dependencies/main.yml | sed '/^$/d' > ./environments/base-p3.10.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.11.yml ./dependencies/binary.yml ./dependencies/main.yml | sed '/^$/d' > ./environments/base-p3.11.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.12.yml ./dependencies/binary.yml ./dependencies/main.yml | sed '/^$/d' > ./environments/base-p3.12.yml


# Viz
./scripts/merge_environment_yaml.py ./dependencies/python3.9.yml  ./dependencies/binary.yml ./dependencies/main.yml ./dependencies/viz.yml | sed '/^$/d' > ./environments/viz-p3.9.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.10.yml ./dependencies/binary.yml ./dependencies/main.yml ./dependencies/viz.yml | sed '/^$/d' > ./environments/viz-p3.10.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.11.yml ./dependencies/binary.yml ./dependencies/main.yml ./dependencies/viz.yml | sed '/^$/d' > ./environments/viz-p3.11.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.12.yml ./dependencies/binary.yml ./dependencies/main.yml ./dependencies/viz.yml | sed '/^$/d' > ./environments/viz-p3.12.yml

# Full
./scripts/merge_environment_yaml.py ./dependencies/python3.9.yml  ./dependencies/binary.yml ./dependencies/dev.yml ./dependencies/docs.yml ./dependencies/main.yml ./dependencies/viz.yml | sed '/^$/d' > ./environments/full-p3.9.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.10.yml ./dependencies/binary.yml ./dependencies/dev.yml ./dependencies/docs.yml ./dependencies/main.yml ./dependencies/viz.yml | sed '/^$/d' > ./environments/full-p3.10.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.11.yml ./dependencies/binary.yml ./dependencies/dev.yml ./dependencies/docs.yml ./dependencies/main.yml ./dependencies/viz.yml | sed '/^$/d' > ./environments/full-p3.11.yml
./scripts/merge_environment_yaml.py ./dependencies/python3.12.yml ./dependencies/binary.yml ./dependencies/dev.yml ./dependencies/docs.yml ./dependencies/main.yml ./dependencies/viz.yml | sed '/^$/d' > ./environments/full-p3.12.yml
