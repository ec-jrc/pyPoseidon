import xarray as xr
import pandas as pd

# to be used for name normalization
LONGITUDE_NAMES = {"longitude", "lon", "x", "Lon", "LONGITUDE", "LON", "X"}
LATITUDE_NAMES = {"latitude", "lat", "y", "Lat", "LATITUDE", "LAT", "Y"}


def normalize_coord_names(dataset: xr.Dataset) -> xr.Dataset:
    """Return a dataset with coords containing "longitude" and "latitude" """
    coords = set(dataset.coords.keys())
    # longitude
    for lon_name in LONGITUDE_NAMES:
        if lon_name in coords:
            break
    else:
        raise ValueError(f"Couldn't normalize longitude: {coords}")
    # latitude
    for lat_name in LATITUDE_NAMES:
        if lat_name in coords:
            break
    else:
        raise ValueError(f"Couldn't normalize latitude: {coords}")
    dataset = dataset.rename({lon_name: "longitude", lat_name: "latitude"})
    return dataset


def normalize_column_names(dataframe: pd.DataFrame) -> pd.DataFrame:
    """Return a dataframe with columns containing "longitude" and "latitude" """
    cols = dataframe.columns.tolist()
    # longitude
    for lon_name in LONGITUDE_NAMES:
        if lon_name in cols:
            break
    else:
        raise ValueError(f"Couldn't normalize longitude: {cols}")
    # latitude
    for lat_name in LATITUDE_NAMES:
        if lat_name in cols:
            break
    else:
        raise ValueError(f"Couldn't normalize latitude: {cols}")

    dic = {lon_name: "longitude", lat_name: "latitude"}
    new_cols = [dic.get(n, n) for n in cols]

    dataframe.columns = new_cols
    return dataframe
