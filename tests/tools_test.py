import os
import pathlib
import shutil

import psutil
import pytest


from pyposeidon import tools


@pytest.mark.skipif(not shutil.which("mpirun"), reason="requires MPI backend")
@pytest.mark.parametrize("use_threads", [True, False])
def test_create_mpirun_script(tmp_path, use_threads):
    target_dir = tmp_path.as_posix()
    cmd = "/bin/schism"
    script_name = "launchSchism.sh"
    script_path = tools.create_mpirun_script(target_dir=target_dir, cmd=cmd, use_threads=use_threads)
    assert script_path == (tmp_path / script_name).as_posix(), "script was not created"
    assert os.access(script_path, os.X_OK), "script is not executable"
    assert not (tmp_path / "outputs").exists(), "outputs subdirectory has not been created"
    content = pathlib.Path(script_path).read_text()
    assert content.endswith(cmd), "the cmd is not being used"
    assert f"-N {psutil.cpu_count(logical=use_threads)}" in content


@pytest.mark.skipif(not shutil.which("mpirun"), reason="requires MPI backend")
@pytest.mark.parametrize("ncores", [1, 2, 44])
def test_create_mpirun_script_ncores(tmp_path, ncores):
    target_dir = tmp_path.as_posix()
    cmd = "/bin/schism"
    script_name = "launchSchism.sh"
    script_path = tools.create_mpirun_script(target_dir=target_dir, cmd=cmd, ncores=ncores)
    assert script_path == (tmp_path / script_name).as_posix(), "script was not created"
    assert os.access(script_path, os.X_OK), "script is not executable"
    assert not (tmp_path / "outputs").exists(), "outputs subdirectory has not been created"
    content = pathlib.Path(script_path).read_text()
    assert content.endswith(cmd), "the cmd is not being used"
    assert f"-N {ncores}" in content
