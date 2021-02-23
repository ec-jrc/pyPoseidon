# Steps for making a new tagged release.

* Stash any uncommitted changes

`git stash`

* bump up. Use minor or patch, e.g.

`bump2version minor`

* Push to origin using the new way

`git push --follow-tags` 

* Go back to the original status

`git stash apply && git stash drop`

* Compile new conda version

`python setup.py clean --all`

if needed...

`export CONDA_BLD_PATH=location_of_conda_channel/`

`conda build conda.recipe --no-test`

`anaconda upload --force location_of_conda_channel/noarch/pyposeidon-X.X.X-py_0.tar.bz2`