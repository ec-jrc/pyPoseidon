import pytest

from pyposeidon import model as pmodel
from pyposeidon.d3d import d3d
from pyposeidon.schism import Schism

from . import DATA_DIR


@pytest.mark.parametrize(
    "model_file,model_class",
    [
        (DATA_DIR / "models" / "d3d_model.json", d3d),
        (DATA_DIR / "models" / "schism_model.json", Schism),
    ],
)
def test_read_from_json(model_file, model_class):
    model = pmodel.read(str(model_file))
    assert isinstance(model, model_class)
