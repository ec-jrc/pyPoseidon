import pandas as pd
import geopandas as gp
import numpy as np
import shapely
from scipy.interpolate import interp1d

# based on https://stackoverflow.com/questions/52014197/how-to-interpolate-a-2d-curve-in-python


def use_spline(df, ds=0.001, method="slinear"):
    points = df.values  # a (nbre_points x nbre_dim) array

    # Linear length along the line:
    distance = np.cumsum(np.sqrt(np.sum(np.diff(points, axis=0) ** 2, axis=1)))

    # Choose n_points
    ips = max(int(distance[-1] / ds), len(distance))

    # normalize
    distance = np.insert(distance, 0, 0) / distance[-1]

    # Interpolation for different methods:
    #    interpolations_methods = ['slinear', 'quadratic', 'cubic']

    alpha = np.linspace(0, 1, ips)

    interpolator = interp1d(distance, points, kind=method, axis=0)
    interpolated_points = interpolator(alpha)

    di = pd.DataFrame(interpolated_points, columns=["lon", "lat"])

    di.iloc[0] = df.iloc[0]
    di.iloc[-1] = df.iloc[-1]

    return di
