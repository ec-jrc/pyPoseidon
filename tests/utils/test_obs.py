from __future__ import annotations

import textwrap

import pandas as pd

from pyposeidon.utils.obs import serialize_stations


def test_serialize_stations(tmp_path):
    expected = textwrap.dedent(
        """\
        1 0 0 0 0 0 0 0 0	 ! https://schism-dev.github.io/schism/master/input-output/optional-inputs.html#stationin-bp-format
        3	 ! number of stations
        1 0.0000000000 0.0000000000 0 	!	 a 0 1.0000000000 1.0000000000 3 157249.3812719440
        2 10.0000000000 5.0000000000 0 	!	 b 1 11.0000000000 4.0000000000 5 157010.1626406018
        3 20.0000000000 0.0000000000 0 	!	 c 2 21.0000000000 1.0000000000 1 157249.3812719441
        """
    )
    stations = pd.DataFrame(
        {
            "lon": [1.0, 11.0, 21.0],
            "lat": [1.0, 4.0, 1.0],
            "unique_id": ["a", "b", "c"],
            "extra_col": ["AA", "BB", "CC"],
            "mesh_index": [0, 1, 2],
            "mesh_lon": [0.0, 10.0, 20.0],
            "mesh_lat": [0.0, 5.0, 0.0],
            "distance": [157249.38127194397, 157010.16264060183, 157249.38127194406],
            "depth": [3, 5, 1],
        }
    )
    path = tmp_path / "station.in"
    serialize_stations(stations, path)
    contents = path.read_text()
    assert contents == expected
