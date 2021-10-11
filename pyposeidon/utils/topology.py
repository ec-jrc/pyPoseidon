import numpy as np


def MakeQuadFaces(Nr, Nc):

    out = np.empty((Nr - 1, Nc - 1, 4), dtype=int)

    r = np.arange(Nr * Nc).reshape(Nr, Nc)

    out[:, :, 0] = r[:-1, :-1]
    out[:, :, 1] = r[:-1, 1:]
    out[:, :, 2] = r[1:, 1:]
    out[:, :, 3] = r[1:, :-1]

    out.shape = (-1, 4)
    return out


def MakeQuadFaces_periodic(Nr, Nc):

    out = np.empty((Nr - 1, Nc, 4), dtype=int)

    r = np.arange(Nr * Nc).reshape(Nr, Nc)

    out[:, :-1, 0] = r[:-1, :-1]
    out[:, :-1, 1] = r[:-1, 1:]
    out[:, :-1, 2] = r[1:, 1:]
    out[:, :-1, 3] = r[1:, :-1]

    out[:, -1, 0] = r[:-1, -1]
    out[:, -1, 1] = r[:-1, 0]
    out[:, -1, 2] = r[1:, 0]
    out[:, -1, 3] = r[1:, -1]

    out.shape = (-1, 4)
    return out


# https://stackoverflow.com/questions/44934631/making-grid-triangular-mesh-quickly-with-numpy
def MakeTriangleFaces(Nr, Nc):

    out = np.empty((Nr - 1, Nc - 1, 2, 3), dtype=int)

    r = np.arange(Nr * Nc).reshape(Nr, Nc)

    out[:, :, 0, 0] = r[:-1, :-1]
    out[:, :, 1, 0] = r[:-1, 1:]
    out[:, :, 0, 1] = r[:-1, 1:]

    out[:, :, 1, 1] = r[1:, 1:]
    out[:, :, :, 2] = r[1:, :-1, None]

    out.shape = (-1, 3)
    return out


def MakeTriangleFaces_periodic(Nr, Nc):

    out = np.empty((Nr - 1, Nc, 2, 3), dtype=int)

    r = np.arange(Nr * Nc).reshape(Nr, Nc)

    out[:, :-1, 0, 0] = r[:-1, :-1]
    out[:, :-1, 1, 0] = r[:-1, 1:]
    out[:, :-1, 0, 1] = r[:-1, 1:]

    out[:, :-1, 1, 1] = r[1:, 1:]
    out[:, :-1, :, 2] = r[1:, :-1, None]

    out[:, -1, 0, 0] = r[:-1, -1]
    out[:, -1, 1, 0] = r[:-1, 0]
    out[:, -1, 0, 1] = r[:-1, 0]
    out[:, -1, 1, 1] = r[1:, 0]
    out[:, -1, :, 2] = r[1:, -1, None]

    out.shape = (-1, 3)
    return out


def quads_to_df(elems, nodes, x="longitude", y="latitude"):
    # cells to polygons
    ap = nodes.loc[elems.a, [x, y]]
    bp = nodes.loc[elems.b, [x, y]]
    cp = nodes.loc[elems.c, [x, y]]
    dp = nodes.loc[elems.d, [x, y]]

    ap["z"] = 0
    bp["z"] = 0
    cp["z"] = 0
    dp["z"] = 0

    elems["ap"] = ap.values.tolist()
    elems["bp"] = bp.values.tolist()
    elems["cp"] = cp.values.tolist()
    elems["dp"] = dp.values.tolist()

    elems["va"] = nodes.loc[elems.a, "d2"].values.tolist()
    elems["vb"] = nodes.loc[elems.b, "d2"].values.tolist()
    elems["vc"] = nodes.loc[elems.c, "d2"].values.tolist()
    elems["vd"] = nodes.loc[elems.d, "d2"].values.tolist()

    return elems.drop(["a", "b", "c", "d"], axis=1)


def quads_to_df_uv(elems, nodes, x="longitude", y="latitude"):
    # cells to polygons
    ap = nodes.loc[elems.a, [x, y]]
    bp = nodes.loc[elems.b, [x, y]]
    cp = nodes.loc[elems.c, [x, y]]
    dp = nodes.loc[elems.d, [x, y]]

    ap["z"] = 0
    bp["z"] = 0
    cp["z"] = 0
    dp["z"] = 0

    elems["ap"] = ap.values.tolist()
    elems["bp"] = bp.values.tolist()
    elems["cp"] = cp.values.tolist()
    elems["dp"] = dp.values.tolist()

    elems["va"] = nodes.loc[elems.a, "d2"].values.tolist()
    elems["vb"] = nodes.loc[elems.b, "d2"].values.tolist()
    elems["vc"] = nodes.loc[elems.c, "d2"].values.tolist()
    elems["vd"] = nodes.loc[elems.d, "d2"].values.tolist()

    return elems.drop(["a", "b", "c", "d"], axis=1)


def quads_to_df_3d(elems, nodes, x="x", y="y", z="z"):
    # cells to polygons
    ap = nodes.loc[elems.a, [x, y, z]]
    bp = nodes.loc[elems.b, [x, y, z]]
    cp = nodes.loc[elems.c, [x, y, z]]
    dp = nodes.loc[elems.d, [x, y, z]]

    elems["ap"] = ap.values.tolist()
    elems["bp"] = bp.values.tolist()
    elems["cp"] = cp.values.tolist()
    elems["dp"] = dp.values.tolist()

    elems["va"] = nodes.loc[elems.a, "d2"].values.tolist()
    elems["vb"] = nodes.loc[elems.b, "d2"].values.tolist()
    elems["vc"] = nodes.loc[elems.c, "d2"].values.tolist()
    elems["vd"] = nodes.loc[elems.d, "d2"].values.tolist()

    return elems.drop(["a", "b", "c", "d"], axis=1)


def tria_to_df_3d(elems, nodes, x="x", y="y", z="z"):
    # cells to polygons
    ap = nodes.loc[elems.a, [x, y, z]]
    bp = nodes.loc[elems.b, [x, y, z]]
    cp = nodes.loc[elems.c, [x, y, z]]

    elems["ap"] = ap.values.tolist()
    elems["bp"] = bp.values.tolist()
    elems["cp"] = cp.values.tolist()

    elems["va"] = nodes.loc[elems.a, "d2"].values.tolist()
    elems["vb"] = nodes.loc[elems.b, "d2"].values.tolist()
    elems["vc"] = nodes.loc[elems.c, "d2"].values.tolist()

    return elems.drop(["a", "b", "c"], axis=1)


def tria_to_df(elems, nodes, x="longitude", y="latitude"):
    # cells to polygons
    ap = nodes.loc[elems.a, [x, y]]
    bp = nodes.loc[elems.b, [x, y]]
    cp = nodes.loc[elems.c, [x, y]]

    ap["z"] = 0
    bp["z"] = 0
    cp["z"] = 0

    elems["ap"] = ap.values.tolist()
    elems["bp"] = bp.values.tolist()
    elems["cp"] = cp.values.tolist()

    elems["va"] = nodes.loc[elems.a, "d2"].values.tolist()
    elems["vb"] = nodes.loc[elems.b, "d2"].values.tolist()
    elems["vc"] = nodes.loc[elems.c, "d2"].values.tolist()

    return elems.drop(["a", "b", "c"], axis=1)
