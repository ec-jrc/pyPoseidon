.PHONY: list
list:
	@LC_ALL=C $(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$'

# Variables
mamba = --mamba
linux = --platform linux-64
mac = --platform osx-64
# Binary
deps_p38_mpich_binary   = --file dependencies/binary_deps.yml --file dependencies/python38.yml --file dependencies/schism59-mpich.yml
deps_p39_mpich_binary   = --file dependencies/binary_deps.yml --file dependencies/python39.yml --file dependencies/schism59-mpich.yml
deps_p38_openmpi_binary = --file dependencies/binary_deps.yml --file dependencies/python38.yml --file dependencies/schism59-openmpi.yml
deps_p39_openmpi_binary = --file dependencies/binary_deps.yml --file dependencies/python39.yml --file dependencies/schism59-openmpi.yml
# No viz
deps_p38_mpich_no_viz   = --file dependencies/python_deps.yml ${deps_p38_mpich_binary}
deps_p38_openmpi_no_viz = --file dependencies/python_deps.yml ${deps_p38_openmpi_binary}
# Viz
deps_p38_mpich_viz   = --file dependencies/viz.yml ${deps_p38_mpich_no_viz}
deps_p38_openmpi_viz = --file dependencies/viz.yml ${deps_p38_openmpi_no_viz}
# Dev
deps_p38_mpich_dev   = --file dependencies/dev.yml ${deps_p38_mpich_viz}
deps_p38_openmpi_dev = --file dependencies/dev.yml ${deps_p38_openmpi_viz}

# Binary Conda locks
conda_lock_linux_p38_mpich_binary:
	conda-lock ${mamba} --check-input-hash ${linux} ${deps_p38_mpich_binary} --filename-template locks/conda-ubuntu-64-p38-mpich-binary.lock
conda_lock_linux_p39_mpich_binary:
	conda-lock ${mamba} --check-input-hash ${linux} ${deps_p39_mpich_binary} --filename-template locks/conda-ubuntu-64-p39-mpich-binary.lock
conda_lock_linux_p38_openmpi_binary:
	conda-lock ${mamba} --check-input-hash ${linux} ${deps_p38_openmpi_binary} --filename-template locks/conda-ubuntu-64-p38-openmpi-binary.lock
conda_lock_linux_p39_openmpi_binary:
	conda-lock ${mamba} --check-input-hash ${linux} ${deps_p39_openmpi_binary} --filename-template locks/conda-ubuntu-64-p39-openmpi-binary.lock
conda_lock_mac_p38_mpich_binary:
	conda-lock ${mamba} --check-input-hash ${mac} ${deps_p38_mpich_binary} --filename-template locks/conda-Macos-64-p38-mpich-binary.lock
conda_lock_mac_p39_mpich_binary:
	conda-lock ${mamba} --check-input-hash ${mac} ${deps_p39_mpich_binary} --filename-template locks/conda-Macos-64-p39-mpich-binary.lock
conda_lock_mac_p38_openmpi_binary:
	conda-lock ${mamba} --check-input-hash ${mac} ${deps_p38_openmpi_binary} --filename-template locks/conda-Macos-64-p38-openmpi-binary.lock
conda_lock_mac_p39_openmpi_binary:
	conda-lock ${mamba} --check-input-hash ${mac} ${deps_p39_openmpi_binary} --filename-template locks/conda-Macos-64-p39-openmpi-binary.lock

# no_viz Conda locks
conda_lock_linux_p38_mpich_no_viz:
	conda-lock ${mamba} --check-input-hash ${linux} ${deps_p38_mpich_no_viz} --filename-template locks/conda-ubuntu-64-p38-mpich-no_viz.lock
conda_lock_linux_p38_openmpi_no_viz:
	conda-lock ${mamba} --check-input-hash ${linux} ${deps_p38_openmpi_no_viz} --filename-template locks/conda-ubuntu-64-p38-openmpi-no_viz.lock
conda_lock_mac_p38_mpich_no_viz:
	conda-lock ${mamba} --check-input-hash ${mac} ${deps_p38_mpich_no_viz} --filename-template locks/conda-Macos-64-p38-mpich-no_viz.lock
conda_lock_mac_p38_openmpi_no_viz:
	conda-lock ${mamba} --check-input-hash ${mac} ${deps_p38_openmpi_no_viz} --filename-template locks/conda-Macos-64-p38-openmpi-no_viz.lock

# viz Conda locks
conda_lock_linux_p38_mpich_viz:
	conda-lock ${mamba} --check-input-hash ${linux} ${deps_p38_mpich_viz} --filename-template locks/conda-ubuntu-64-p38-mpich-viz.lock
conda_lock_linux_p38_openmpi_viz:
	conda-lock ${mamba} --check-input-hash ${linux} ${deps_p38_openmpi_viz} --filename-template locks/conda-ubuntu-64-p38-openmpi-viz.lock
conda_lock_mac_p38_mpich_viz:
	conda-lock ${mamba} --check-input-hash ${mac} ${deps_p38_mpich_viz} --filename-template locks/conda-Macos-64-p38-mpich-viz.lock
conda_lock_mac_p38_openmpi_viz:
	conda-lock ${mamba} --check-input-hash ${mac} ${deps_p38_openmpi_viz} --filename-template locks/conda-Macos-64-p38-openmpi-viz.lock

# Dev Conda locks
conda_lock_linux_p38_mpich_dev:
	conda-lock ${mamba} --check-input-hash ${linux} ${deps_p38_mpich_dev} --filename-template locks/conda-ubuntu-64-p38-mpich-dev.lock
conda_lock_linux_p38_openmpi_dev:
	conda-lock ${mamba} --check-input-hash ${linux} ${deps_p38_openmpi_dev} --filename-template locks/conda-ubuntu-64-p38-openmpi-dev.lock
conda_lock_mac_p38_mpich_dev:
	conda-lock ${mamba} --check-input-hash ${mac} ${deps_p38_mpich_dev} --filename-template locks/conda-Macos-64-p38-mpich-dev.lock
conda_lock_mac_p38_openmpi_dev:
	conda-lock ${mamba} --check-input-hash ${mac} ${deps_p38_openmpi_dev} --filename-template locks/conda-Macos-64-p38-openmpi-dev.lock

conda_lock_linux: \
	conda_lock_linux_p38_mpich_binary \
	conda_lock_linux_p39_mpich_binary \
	conda_lock_linux_p38_openmpi_binary \
	conda_lock_linux_p39_openmpi_binary \
	conda_lock_linux_p38_mpich_no_viz \
	conda_lock_linux_p38_mpich_viz \
	conda_lock_linux_p38_mpich_dev \
	conda_lock_linux_p38_openmpi_no_viz \
	conda_lock_linux_p38_openmpi_viz \
	conda_lock_linux_p38_openmpi_dev

conda_lock_mac: \
	conda_lock_mac_p38_mpich_binary \
	conda_lock_mac_p39_mpich_binary \
	conda_lock_mac_p38_openmpi_binary \
	conda_lock_mac_p39_openmpi_binary \
	conda_lock_mac_p38_mpich_no_viz \
	conda_lock_mac_p38_mpich_viz \
	conda_lock_mac_p38_mpich_dev \
	conda_lock_mac_p38_openmpi_no_viz \
	conda_lock_mac_p38_openmpi_viz \
	conda_lock_mac_p38_openmpi_dev

conda_lock: \
	conda_lock_linux \
	conda_lock_mac

poetry_lock:
	poetry lock --no-update
	poetry export --without-hashes -f requirements.txt -o locks/requirements.txt
	poetry export --without-hashes -f requirements.txt --extras viz -o locks/requirements-viz.txt
	poetry export --without-hashes -f requirements.txt --extras testing -o locks/requirements-testing.txt
	poetry export --without-hashes -f requirements.txt --extras testing --extras viz --extras docs --dev -o locks/requirements-full.txt

lock: \
	poetry_lock \
	conda_lock
