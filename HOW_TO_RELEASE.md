# Steps for making a new tagged release.

* Stash any uncommitted changes

`git stash`

* bump up. Use minor or patch, e.g.

`bump2version minor/patch`

- check for issues

`git show`

- After fixing issues

`git add .`

`git commit --amend`

## Build pip version

`python3 -m venv .venv`

`source .venv/bin/activate`

`poetry update`

- clean dist/

`poetry build`


- test CI on github

`git push`

## Compile anaconda conda version

`python setup.py clean --all`

- if needed `export CONDA_BLD_PATH=location_of_conda_channel/`

`conda build conda.recipe`

`anaconda upload --force location_of_conda_channel/noarch/pyposeidon-X.X.X-py_0.tar.bz2`

## Tag

* If all works create tag and push tag e.g.

`git tag -a -m 'Release: 0.4.2 -> 0.5.0' '0.5.0'`

`git push --follow-tags`


## Set release on Github


## Update conda-forge.

- Clone https://github.com/conda-forge/pyposeidon-feedstock and update the version number and sha256 in meta.yaml (On OS X, you can calculate sha256 with `shasum -a 256 file.tar.gz`).

- Submit a pull request (and merge it, once CI passes).


## Publish on pypi

`poetry publish`

## Go back to the original status

`git stash apply && git stash drop`
