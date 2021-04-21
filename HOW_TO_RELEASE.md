# Steps for making a new tagged release.

* Stash any uncommitted changes

`git stash`

* bump up. Use minor or patch, e.g.

`bump2version minor/patch`

* Push to origin using the new way

`git push --follow-tags`


## Build pip version

`python3 -m venv .venv`

`source .venv/bin/activate`

`poetry install`

- check in pyproject.toml for issues

`poetry update`

- clean dist/

`poetry build`

`poetry publish`

## Update conda-forge.

- Clone https://github.com/conda-forge/pyposeidon-feedstock and update the version number and sha256 in meta.yaml (On OS X, you can calculate sha256 with `shasum -a 256 file.tar.gz`).

- Submit a pull request (and merge it, once CI passes).


## Compile anaconda conda version

`python setup.py clean --all`

- if needed `export CONDA_BLD_PATH=location_of_conda_channel/`

`conda build conda.recipe --no-test`

`anaconda upload --force location_of_conda_channel/noarch/pyposeidon-X.X.X-py_0.tar.bz2`


## Go back to the original status

`git stash apply && git stash drop`
