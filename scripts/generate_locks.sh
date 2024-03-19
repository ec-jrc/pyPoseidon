#!/usr/bin/env bash
#

set -xeuo pipefail

# Schism
conda-lock lock --mamba --check-input-hash -p linux-64 -p osx-64 -f dependencies/schism_openmpi.yml 	--lockfile locks/schism_openmpi.yml
conda-lock lock --mamba --check-input-hash -p linux-64 -p osx-64 -f dependencies/schism_mpich.yml 		--lockfile locks/schism_mpich.yml
conda-lock render -p linux-64 --filename-template locks/conda-ubuntu-latest-schism_openmpi.lock locks/schism_openmpi.yml
conda-lock render -p linux-64 --filename-template locks/conda-ubuntu-latest-schism_mpich.lock   locks/schism_mpich.yml
conda-lock render -p osx-64 --filename-template locks/conda-macos-latest-schism_openmpi.lock    locks/schism_openmpi.yml
conda-lock render -p osx-64 --filename-template locks/conda-macos-latest-schism_mpich.lock      locks/schism_mpich.yml

# delft3d
conda-lock lock --mamba --check-input-hash -p linux-64 -p osx-64 -f dependencies/delft3d_openmpi.yml  --lockfile locks/delft3d_openmpi.yml
conda-lock lock --mamba --check-input-hash -p linux-64 -p osx-64 -f dependencies/delft3d_mpich.yml    --lockfile locks/delft3d_mpich.yml
conda-lock render -p linux-64 --filename-template locks/conda-ubuntu-latest-delft3d_openmpi.lock  locks/delft3d_openmpi.yml
conda-lock render -p linux-64 --filename-template locks/conda-ubuntu-latest-delft3d_mpich.lock    locks/delft3d_mpich.yml
conda-lock render -p osx-64 --filename-template locks/conda-macos-latest-delft3d_openmpi.lock     locks/delft3d_openmpi.yml
conda-lock render -p osx-64 --filename-template locks/conda-macos-latest-delft3d_mpich.lock       locks/delft3d_mpich.yml

# binary
conda-lock lock --mamba --check-input-hash -p linux-64 -p osx-64 -f environments/binary-p3.9.yml    --lockfile locks/binary-p3.9.yml
conda-lock lock --mamba --check-input-hash -p linux-64 -p osx-64 -f environments/binary-p3.10.yml   --lockfile locks/binary-p3.10.yml
conda-lock lock --mamba --check-input-hash -p linux-64 -p osx-64 -f environments/binary-p3.11.yml   --lockfile locks/binary-p3.11.yml
conda-lock lock --mamba --check-input-hash -p linux-64 -p osx-64 -f environments/binary-p3.12.yml   --lockfile locks/binary-p3.12.yml
conda-lock render -p linux-64 --filename-template locks/conda-ubuntu-latest-binary-p3.9.lock    locks/binary-p3.9.yml
conda-lock render -p linux-64 --filename-template locks/conda-ubuntu-latest-binary-p3.10.lock   locks/binary-p3.10.yml
conda-lock render -p linux-64 --filename-template locks/conda-ubuntu-latest-binary-p3.11.lock   locks/binary-p3.11.yml
conda-lock render -p linux-64 --filename-template locks/conda-ubuntu-latest-binary-p3.12.lock   locks/binary-p3.12.yml
conda-lock render -p osx --filename-template locks/conda-macos-latest-binary-p3.9.lock    locks/binary-p3.9.yml
conda-lock render -p osx --filename-template locks/conda-macos-latest-binary-p3.10.lock   locks/binary-p3.10.yml
conda-lock render -p osx --filename-template locks/conda-macos-latest-binary-p3.11.lock   locks/binary-p3.11.yml
conda-lock render -p osx --filename-template locks/conda-macos-latest-binary-p3.12.lock   locks/binary-p3.12.yml

# Telemac
conda-lock lock --mamba --check-input-hash -p linux-64 -f environments/binary-telemac-p3.9.yml 	  --lockfile locks/binary-telemac-p3.9.yml
conda-lock lock --mamba --check-input-hash -p linux-64 -f environments/binary-telemac-p3.10.yml 	--lockfile locks/binary-telemac-p3.10.yml
conda-lock lock --mamba --check-input-hash -p linux-64 -f environments/binary-telemac-p3.11.yml 	--lockfile locks/binary-telemac-p3.11.yml
conda-lock lock --mamba --check-input-hash -p linux-64 -f environments/binary-telemac-p3.12.yml 	--lockfile locks/binary-telemac-p3.12.yml
conda-lock render -p linux-64 --filename-template locks/conda-ubuntu-latest-binary-telemac-p3.9.lock   locks/binary-telemac-p3.9.yml
conda-lock render -p linux-64 --filename-template locks/conda-ubuntu-latest-binary-telemac-p3.10.lock  locks/binary-telemac-p3.10.yml
conda-lock render -p linux-64 --filename-template locks/conda-ubuntu-latest-binary-telemac-p3.11.lock  locks/binary-telemac-p3.11.yml
conda-lock render -p linux-64 --filename-template locks/conda-ubuntu-latest-binary-telemac-p3.12.lock  locks/binary-telemac-p3.12.yml


# Pyposeidon-Base
# Pyposeidon-Viz
# Pyposeidon-Full
# conda-lock lock --mamba --check-input-hash -p linux-64 -p osx-64 -f environments/full-p3.9.yml --lockfile locks/full-p3.9.yml
# conda-lock render -p linux-64 --filename-template locks/conda-ubuntu-latest-full-p3.9.lock locks/full-p3.9.yml
# conda-lock render -p osx-64   --filename-template locks/conda-macos-latest-full-p3.9.lock 	locks/full-p3.9.yml
