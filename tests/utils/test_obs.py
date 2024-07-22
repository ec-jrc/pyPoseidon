from __future__ import annotations

import textwrap
import pandas as pd

from pyposeidon.utils.obs import serialize_stations


def test_serialize_stations(tmp_path):
    expected = textwrap.dedent(
        """\
        1 0 0 0 0 0 0 0 0	 ! https://schism-dev.github.io/schism/master/input-output/optional-inputs.html#stationin-bp-format
        3	 ! number of stations
        1 0.0000000000 0.0000000000 0 	!	 0 1.0000000000 1.0000000000 157249.3812719440 AA a
        2 10.0000000000 5.0000000000 0 	!	 1 11.0000000000 4.0000000000 157010.1626406018 BB b
        3 20.0000000000 0.0000000000 0 	!	 2 21.0000000000 1.0000000000 157249.3812719441 CC c
        """
    )
    stations = pd.DataFrame({
        'lon': [1., 11., 21.],
        'lat': [1., 4., 1.],
        'unique_id': ["a", "b", "c"],
        'extra_col': ["AA", "BB", "CC"],
        'mesh_index': [0, 1, 2],
        'mesh_lon': [0., 10., 20.],
        'mesh_lat': [0., 5., 0.],
        'distance': [157249.38127194397, 157010.16264060183, 157249.38127194406],
    })
    path = tmp_path / "station.in"
    serialize_stations(stations, path)
    contents = path.read_text()
    assert contents == expected
