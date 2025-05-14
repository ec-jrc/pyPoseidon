import xarray as xr
import pandas as pd

# to be used for name normalization
LONGITUDE_NAMES = {"longitude", "lon", "x", "Lon", "LONGITUDE", "LON", "X", "SCHISM_hgrid_node_x"}
LATITUDE_NAMES = {"latitude", "lat", "y", "Lat", "LATITUDE", "LAT", "Y", "SCHISM_hgrid_node_y"}
LOCATION_NAMES = {"location", "loc", "LOCATION", "Station_Name", "Station"}

TEL_MAP_LONG = {
    "VELOCITY U      ": "u-current-velocity",
    "VELOCITY V      ": "v-current-velocity",
    "FREE SURFACE    ": "free-surface",
    "WAVE HEIGHT HM0 ": "significant-wave-height",
    "MEAN DIRECTION  ": "mean-direction",
    "MEAN PERIOD TM01": "mean-period-tm01",
    "MEAN PERIOD TM02": "mean-period-tm02",
    "MEAN PERIOD TMOY": "mean-period-tmoy",
    "PEAK PERIOD TPR8": "peak-period-tpr8",
    "PEAK PERIOD TPD ": "peak-period-tpd",
    "WIND VELOCITY U ": "10m-u-component-of-wind",
    "WIND VELOCITY V ": "10m-v-component-of-wind",
    "SURFACE PRESSURE": "mean-sea-level-pressure",
    "AIR TEMPERATURE ": "air-temperature",
}
TEL_MAP_SHORT = {
    "VELOCITY U      ": "u",
    "VELOCITY V      ": "v",
    "FREE SURFACE    ": "elev",
    "WAVE HEIGHT HM0 ": "hm0",
    "MEAN DIRECTION  ": "dir",
    "MEAN PERIOD TM01": "tm01",
    "MEAN PERIOD TM02": "tm02",
    "MEAN PERIOD TMOY": "tmoy",
    "PEAK PERIOD TPR8": "tp",
    "PEAK PERIOD TPD ": "tp",
    "WIND VELOCITY U ": "u10",
    "WIND VELOCITY V ": "v10",
    "SURFACE PRESSURE": "msl",
    "AIR TEMPERATURE ": "tair",
}

TEL_MAP_SHORT2 = {
    "x": "longitude",
    "y": "latitude",
    "U": "u",
    "V": "v",
    "S": "elev",
    "Z": "elev",
    "WH": "hm0",
    "DIRM": "dir",
    "TM01": "tm01",
    "TM02": "tm02",
    "TMOY": "tmoy",
    "TPR8": "tp",
    "TPD": "tp",
    "WINDX": "u10",
    "WINDY": "v10",
    "PATM": "msl",
    "TAIR": "tair",
}


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
    # location
    for location in LOCATION_NAMES:
        if location in cols:
            break
    else:
        raise ValueError(f"Couldn't normalize location: {cols}")
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

    dic = {lon_name: "longitude", lat_name: "latitude", location: "location"}
    new_cols = [dic.get(n, n) for n in cols]

    dataframe.columns = new_cols
    return dataframe


def normalize_varnames(varnames, mapping=TEL_MAP_SHORT2, mapping1=TEL_MAP_SHORT):
    normalized_names = []
    for v in varnames:
        if v in mapping:
            normalized_names.append(mapping[v])
        elif v in mapping1:
            normalized_names.append(mapping1[v])
        else:
            normalized_names.append(v.strip().lower())
    return normalized_names
