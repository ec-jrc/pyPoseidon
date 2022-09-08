import pandas as pd
import xarray as xr

from pyposeidon.utils.topology import tria_to_df, tria_to_df_3d
from pyposeidon.utils.stereo import stereo_to_3d


def to_sq(df, fpos):

    with open(fpos, "w") as f:
        f.write("//*********************************************************************\n")
        f.write("// *\n")
        f.write("// *  pyposeidon\n")
        f.write("// *\n")
        f.write("// *  Scalar 2D post-processing view\n")
        f.write("// *\n")
        f.write("// *********************************************************************/\n\n")

        f.write("// This view contains a scalar field defined on quads.\n")
        f.write("\n")
        f.write('View "{}" {}\n'.format("bgmesh", "{"))
        for idx, vals in df.iterrows():
            f.write(
                "SQ({},{},{},{},{},{},{},{},{},{},{},{}){{{},{},{},{}}};\n".format(
                    vals.ap[0],
                    vals.ap[1],
                    vals.ap[2],
                    vals.bp[0],
                    vals.bp[1],
                    vals.bp[2],
                    vals.cp[0],
                    vals.cp[1],
                    vals.cp[2],
                    vals.dp[0],
                    vals.dp[1],
                    vals.dp[2],
                    vals.va,
                    vals.vb,
                    vals.vc,
                    vals.vd,
                )
            )

        f.write("{};\n".format("}"))


def to_st(df, fpos):

    with open(fpos, "w") as f:
        f.write("//*********************************************************************\n")
        f.write("// *\n")
        f.write("// *  pyposeidon\n")
        f.write("// *\n")
        f.write("// *  Scalar 2D post-processing view\n")
        f.write("// *\n")
        f.write("// *********************************************************************/\n\n")

        f.write("// This view contains a scalar field defined on triangles.\n")
        f.write("\n")
        f.write('View "{}" {}\n'.format("bgmesh", "{"))
        for idx, vals in df.iterrows():
            f.write(
                "ST({},{},{},{},{},{},{},{},{}){{{},{},{}}};\n".format(
                    vals.ap[0],
                    vals.ap[1],
                    vals.ap[2],
                    vals.bp[0],
                    vals.bp[1],
                    vals.bp[2],
                    vals.cp[0],
                    vals.cp[1],
                    vals.cp[2],
                    vals.va,
                    vals.vb,
                    vals.vc,
                )
            )

        f.write("{};\n".format("}"))


def to_global_pos(nodes, elems, fpos, **kwargs):

    use_bindings = kwargs.get("use_bindings", True)
    R = kwargs.get("R", 1.0)

    if use_bindings:
        # Keep stereographic

        elems["d"] = 0

        sv = 4 * R**2 / (nodes.u**2 + nodes.v**2 + 4 * R**2)
        nodes["d2"] = nodes.d2 / sv

        dout = tria_to_df(elems, nodes, x="u", y="v")
        # save bgmesh
        to_st(dout, fpos)

    else:
        # use 3D
        x, y, z = stereo_to_3d(nodes.u.values, nodes.v.values)

        nodes["x"] = x
        nodes["y"] = y
        nodes["z"] = z

        # create output dataframe
        dout = tria_to_df_3d(elems, nodes)
        # save bgmesh
        to_st(dout, fpos)

    ## make dataset
    els = xr.DataArray(
        elems.loc[:, ["a", "b", "c"]],
        dims=["nSCHISM_hgrid_face", "nMaxSCHISM_hgrid_face_nodes"],
        name="SCHISM_hgrid_face_nodes",
    )

    nod = (
        nodes.loc[:, ["u", "v"]]
        .to_xarray()
        .rename(
            {
                "index": "nSCHISM_hgrid_node",
                "u": "SCHISM_hgrid_node_x",
                "v": "SCHISM_hgrid_node_y",
            }
        )
    )
    nod = nod.drop_vars("nSCHISM_hgrid_node")

    bg = xr.Dataset({"h": (["nSCHISM_hgrid_node"], nodes.d2.values)})

    dh = xr.merge([nod, els, bg])

    return dh
