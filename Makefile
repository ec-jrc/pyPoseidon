.PHONY: list
list:
	@LC_ALL=C $(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | grep -E -v -e '^[^[:alnum:]]' -e '^$@$$'

conda_lock: \
	generate_envs
	conda-lock lock --mamba --check-input-hash -p linux-64 -p osx-64 -f environments/binary_mpich.yml --lockfile locks/binary_mpich.yml
	conda-lock lock --mamba --check-input-hash -p linux-64 -p osx-64 -f environments/binary_openmpi.yml --lockfile locks/binary_openmpi.yml
	conda-lock lock --mamba --check-input-hash -p linux-64 -p osx-64 -f environments/base.yml --lockfile locks/base.yml
	conda-lock lock --mamba --check-input-hash -p linux-64 -p osx-64 -f environments/viz.yml --lockfile locks/viz.yml
	conda-lock lock --mamba --check-input-hash -p linux-64 -p osx-64 -f environments/full.yml --lockfile locks/full.yml
	#
	conda-lock render -p linux-64 --filename-template locks/conda-ubuntu-64-p38-openmpi-binary.lock 			locks/binary_openmpi.yml
	conda-lock render -p linux-64 --filename-template locks/conda-ubuntu-64-p38-mpich-binary.lock   			locks/binary_mpich.yml
	conda-lock render -p linux-64 --filename-template locks/conda-ubuntu-64-p38-mpich-base.lock           locks/base.yml
	conda-lock render -p linux-64 --filename-template locks/conda-ubuntu-64-p38-mpich-viz.lock            locks/viz.yml
	conda-lock render -p linux-64 --filename-template locks/conda-ubuntu-64-p38-mpich-full.lock           locks/full.yml
	#
	conda-lock render -p osx-64 --filename-template locks/conda-macos-64-p38-openmpi-binary.lock  			locks/binary_openmpi.yml
	conda-lock render -p osx-64 --filename-template locks/conda-macos-64-p38-mpich-binary.lock    			locks/binary_mpich.yml
	conda-lock render -p osx-64 --filename-template locks/conda-macos-64-p38-mpich-base.lock            locks/base.yml
	conda-lock render -p osx-64 --filename-template locks/conda-macos-64-p38-mpich-viz.lock             locks/viz.yml
	conda-lock render -p osx-64 --filename-template locks/conda-macos-64-p38-mpich-full.lock            locks/full.yml
	conda-lock render -p osx-64 --filename-template locks/conda-macos-64-p38-mpich-full.lock            locks/full.yml

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
	poetry export --without-hashes -f requirements.txt --extras testing --extras viz --extras docs --with dev -o locks/requirements-full.txt

lock: \
	poetry_lock \
	conda_lock
