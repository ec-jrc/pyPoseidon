"""Get boundaries"""
import numpy as np
from pyposeidon.tide import *


def grouper(iterable):
    prev = None
    group = []
    for item in iterable:
        if not prev or item - prev <= 1:
            group.append(item)
        else:
            yield group
            group = [item]
        prev = item
    if group:
        yield group


class Bound:
    def __init__(self, **kwargs):
        btype = kwargs.get("btype", None)


class Box(Bound):
    def __init__(self, **kwargs):

        lons = kwargs.get("lons", None)
        lats = kwargs.get("lats", None)
        ba = kwargs.get("dem", None)
        n = kwargs.get("cn", None)
        tpath = kwargs.get("tpath", "./")
        tmodel = kwargs.get("tmodel", "./")

        dic = {
            "West": ba[:, 0],
            "East": ba[:, -1],
            "South": ba[0, :],
            "North": ba[-1, :],
        }
        idic = {
            "West": [0, -99],
            "East": [lons.shape[0], -99],
            "South": [-99, 0],
            "North": [-99, lats.shape[0]],
        }

        for bound in ["West", "East", "North", "South"]:
            chunks = []

            bv = dic[bound]  # Identify boundary

            v = np.isfinite(bv).nonzero()  # get non zero bathymetry part
            branches = dict(enumerate(grouper(list(v[0])), 1))

            # iterate over all branches
            s = []
            for b in branches:
                if len(branches[b]) > 0:
                    s = branches[b][:] if branches[b][0] not in [0, 1] else branches[b][1:]

            # make chunks of n points if range is large
            chunks = []
            if len(s) > 0:
                chunks = [s[i : i + n] for i in xrange(0, len(s), n if n < len(s) else len(s))]

            q = []

            for ch in chunks:

                k1 = [ch[0] if x == -99 else x for x in idic[bound]]
                k2 = [ch[-1] if x == -99 else x for x in idic[bound]]

                q.append([[k1[0], k1[1]], [k2[0], k2[1]]])

            # store the values
            setattr(self, bound, q)
