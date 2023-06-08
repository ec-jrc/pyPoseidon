def get_dem_contours(blevels, dem):
    if len(blevels) == 1:
        dvalues = [blevels[0] - 1, blevels[0]]
    elif len(blevels) == 2:
        dvalues = blevels
    else:
        logger.error("blevels not set properly")
        sys.exit(1)

    # check if in range
    r_check = [dem.elevation.values.min() < x < dem.elevation.values.max() for x in dvalues]

    if not any(r_check):
        logger.warning(
            "blevels not within DEM range {}:{}".format(dem.elevation.values.min(), dem.elevation.values.max())
        )
        return gp.GeoDataFrame(geometry=[])
    try:
        dem.adjusted.plot.contour(levels=dvalues)
    except:
        dem.elevation.plot.contour(levels=dvalues)

    # get graph
    cs = plt.gca()

    # outer contour
    outc = []
    outl = cs.collections[0].get_paths()
    for l in range(len(outl)):
        v = outl[l].vertices
        try:
            outc.append(shapely.geometry.LineString(v))
        except:
            pass

    # inner contour
    inc = []
    inl = cs.collections[1].get_paths()
    for l in range(len(inl)):
        v = inl[l].vertices
        try:
            inc.append(shapely.geometry.LineString(v))
        except:
            pass

    # Get continuous area
    lon_min = dem.longitude.min()
    lon_max = dem.longitude.max()
    lat_min = dem.latitude.min()
    lat_max = dem.latitude.max()

    # create a polygon of the lat/lon window
    grp = shapely.geometry.Polygon([(lon_min, lat_min), (lon_min, lat_max), (lon_max, lat_max), (lon_max, lat_min)])

    # outer contour
    b1 = []
    for geo_ in outc:
        ee = shapely.ops.split(grp, geo_)
        try:
            b1.append(ee[1])
        except:
            pass
    # inner contour
    b2 = []
    for geo_ in inc:
        ee = shapely.ops.split(grp, geo_)
        try:
            b2.append(ee[1])
        except:
            pass

    w1 = gp.GeoDataFrame(geometry=b1)
    w1_ = w1.unary_union
    w1 = gp.GeoDataFrame(geometry=[w1_])

    w2 = gp.GeoDataFrame(geometry=b2)
    w2_ = w2.unary_union
    w2 = gp.GeoDataFrame(geometry=[w2_])

    plt.close("all")
    return w1.loc[~w1.is_empty]  # cleanup
