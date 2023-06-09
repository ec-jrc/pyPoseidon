.PHONY: list
list:
	@LC_ALL=C $(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | grep -E -v -e '^[^[:alnum:]]' -e '^$@$$'

init:
	poetry install -E viz --with dev --with docs --with test

test:
	pytest -vl --durations=10 -m 'not (schism or delft or viz or slow)' -n auto

test_fail:
	pytest -vl --lf --runschism --rundelft --runviz

test_viz:
	pytest -vl --durations=10 --runviz -m viz

test_schism:
	pytest -vl --durations=10 --runschism -m schism

test_delft:
	pytest -vl --durations=10 --rundelft -m delft

test_full:
	pytest -vl --durations=20 --runschism --rundelft --runviz

conda_lock:
	./scripts/generate_envs.sh
	./scripts/generate_locks.sh

poetry_lock:
	poetry lock --no-update
	poetry export --without-hashes -f requirements.txt -o locks/requirements.txt
	poetry export --without-hashes -f requirements.txt --extras viz -o locks/requirements-viz.txt
	poetry export --without-hashes -f requirements.txt --with dev --with test -o locks/requirements-ci.txt
	poetry export --without-hashes -f requirements.txt --extras viz --with dev --with docs --with test -o locks/requirements-full.txt

lock: \
	poetry_lock \
	conda_lock
