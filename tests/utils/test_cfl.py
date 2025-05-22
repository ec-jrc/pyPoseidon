from pyposeidon.utils.cfl import parse_gr3

from .. import DATA_DIR


def test_parse_gr3_with_boundaries():
    parsed = parse_gr3(DATA_DIR / "hgrid_with_multiple_boundaries.gr3", include_boundaries=True)
    assert "nodes" in parsed
    assert len(parsed["nodes"]) == 1165
    assert "elements" in parsed
    assert len(parsed["elements"]) == 1949
    assert "boundaries" in parsed
    assert len(parsed["boundaries"]) == 2
    assert "open" in parsed["boundaries"]
    assert len(parsed["boundaries"]["open"]) == 2
    assert 1 in parsed["boundaries"]
    assert len(parsed["boundaries"][1]) == 2
