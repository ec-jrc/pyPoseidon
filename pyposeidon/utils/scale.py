import numpy as np
import xarray as xr


def scale_dem(b, res_min, res_max, **kwargs):

    #    b.columns = ["z"]

    b.loc[b.z >= -10, "z"] = -1.0e-4  # normalize to only negative values

    b.z = np.sqrt(-b.z) / 0.5  # scale

    # adjust scale

    bg = b.z.values

    a2 = (bg - bg.min()) / (bg.max() - bg.min())

    d2 = res_min + a2 * (res_max - res_min)

    b["d2"] = d2

    nodes = b.reset_index(drop=True)

    nodes["z"] = 0

    #    subspace = kwargs.get('subspace', None)

    #    if subspace is not None:
    #        mask = ...
    #        hfun[mask] =

    return nodes
