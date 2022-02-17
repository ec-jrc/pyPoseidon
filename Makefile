.PHONY: list
list:
	@LC_ALL=C $(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$'


# Variables
conda_lock_cmd = conda-lock --mamba --check-input-hash
linux = --platform linux-64
macos = --platform osx-64

conda_lock_linux:
	${conda_lock_cmd} ${linux} --file environments/binary_mpich.yml --filename-template locks/conda-ubuntu-64-p38-mpich-binary.lock
	${conda_lock_cmd} ${linux} --file environments/binary_openmpi.yml --filename-template locks/conda-ubuntu-64-p38-openmpi-binary.lock
	${conda_lock_cmd} ${linux} --file environments/base.yml --filename-template locks/conda-ubuntu-64-p38-mpich-base.lock
	${conda_lock_cmd} ${linux} --file environments/viz.yml --filename-template locks/conda-ubuntu-64-p38-mpich-viz.lock
	${conda_lock_cmd} ${linux} --file environments/full.yml --filename-template locks/conda-ubuntu-64-p38-mpich-full.lock

conda_lock_macos:
	${conda_lock_cmd} ${macos} --file environments/binary_mpich.yml --filename-template locks/conda-macos-64-p38-mpich-binary.lock
	${conda_lock_cmd} ${macos} --file environments/binary_openmpi.yml --filename-template locks/conda-macos-64-p38-openmpi-binary.lock
	${conda_lock_cmd} ${macos} --file environments/base.yml --filename-template locks/conda-macos-64-p38-mpich-base.lock
	${conda_lock_cmd} ${macos} --file environments/viz.yml --filename-template locks/conda-macos-64-p38-mpich-viz.lock
	${conda_lock_cmd} ${macos} --file environments/full.yml --filename-template locks/conda-macos-64-p38-mpich-full.lock

conda_lock: \
	generate_envs \
	conda_lock_linux \
	conda_lock_macos

generate_binary_env:
	./scripts/merge_environment_yaml.py \
		./dependencies/binary.yml \
		./dependencies/solvers_mpich.yml \
	> ./environments/binary_mpich.yml
	./scripts/merge_environment_yaml.py \
		./dependencies/binary.yml \
		./dependencies/solvers_openmpi.yml \
	> ./environments/binary_openmpi.yml

generate_base_env:
	./scripts/merge_environment_yaml.py \
		./dependencies/binary.yml \
		./dependencies/python.yml \
		./dependencies/solvers_mpich.yml \
	> ./environments/base.yml

generate_viz_env:
	./scripts/merge_environment_yaml.py \
		./dependencies/binary.yml \
		./dependencies/python.yml \
		./dependencies/solvers_mpich.yml \
		./dependencies/viz.yml \
	> ./environments/viz.yml

generate_full_env:
	scripts/merge_environment_yaml.py \
		./dependencies/binary.yml \
		./dependencies/dev.yml \
		./dependencies/docs.yml \
		./dependencies/python.yml \
		./dependencies/solvers_mpich.yml \
		./dependencies/viz.yml \
	> ./environments/full.yml

generate_envs: \
	generate_binary_env \
	generate_base_env \
	generate_viz_env \
	generate_full_env

poetry_lock:
	poetry lock --no-update
	poetry export --without-hashes -f requirements.txt -o locks/requirements.txt
	poetry export --without-hashes -f requirements.txt --extras viz -o locks/requirements-viz.txt
	poetry export --without-hashes -f requirements.txt --extras testing --extras viz --extras docs --dev -o locks/requirements-full.txt

lock: \
	poetry_lock \
	conda_lock
