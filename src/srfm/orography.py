"""Works with Earth's orography.

Plots orography, calculates elevation at given location.

- Name: orography
- Parent package: srfm
- Author: Antonin Knizek
- Contributors:
- Date: 5 August 2025
"""

import netCDF4
import matplotlib.pyplot as plt
import numpy as np
from cartopy import crs as ccrs, feature as cfeature
from scipy.interpolate import RegularGridInterpolator
from cartopy.mpl.ticker import (
    LatitudeFormatter,
    LongitudeFormatter,
    LatitudeLocator,
    LongitudeLocator,
)
import matplotlib.style as mplstyle
from importlib.resources import files, as_file

"""
This is the default file containing orographical data.
To use a different file, change here.
Mind that the code was written to work with this file only.
"""


def load_orography():
    """Loads orography.nc file.

    The file is part of SRFM.

    """
    with as_file(files("srfm.data") / "orography.nc") as path:
        ds = netCDF4.Dataset(path, "r")
    return ds


# data = netCDF4.Dataset("/home/k/knizek/Documents/srfm_main/srfm/data/orography.nc", "r")


def plot_Earth():
    """Plots the Earth's orography on a map."""
    data = load_orography()

    x, y = np.meshgrid(data.variables["lon"][:], data.variables["lat"][:])

    plt.figure(figsize=(11, 8.5))
    projPC = ccrs.PlateCarree()
    ax = plt.subplot(1, 1, 1, projection=projPC)

    ax.coastlines()

    ##plor cloud cover
    cmap = "Grays"
    ax.scatter(
        x, y, c=data.variables["orog"][:], transform=ccrs.PlateCarree(), cmap=cmap
    )
    # set ticks automatically:
    gl = ax.gridlines(
        draw_labels=["bottom", "left"],  # draws labels
        dms=False,
    )

    # sets tick spacing:
    # gl.xlocator = mticker.MultipleLocator(1)
    #    gl.xlocator = LongitudeLocator(nbins=180) #for some reason sets the spacing as 360/nbins, if not specified, creates ticks every 36 deg, which is useless for zoomed plots
    #    gl.ylocator = LatitudeLocator()
    gl.xlines = False
    gl.ylines = False

    ax.set_xlabel("Longitude")
    ax.set_ylabel("Lattitude")

    mplstyle.use("fast")
    plt.show()
    return


def get_elevation(coords):
    """Calculates elevation at given coordinates.

    Args:
        coords (array-like): Lat, lon pairs at which to evaluate elevation.
                             2D array expected.
                             For single point, 1D array or list [lat,lon] are accepted
                             and converted to 2D array.
                             First column ([:,0]) is latitudes,
                             second column ([:,1]) is longitudes.
                             Lat must be within (-89.516,,89.866).
                             Lon must be within (-180,179.473).

    Returns:
        ele (array-like): Elevation array. Shape coords.shape[0]. Units [m].

    """

    data = load_orography()

    if isinstance(coords, list) and len(coords) == 2:
        coords = np.array(coords)

    if isinstance(coords, np.ndarray) and coords.ndim == 1:
        coords = np.expand_dims(coords, axis=0)

    interp = RegularGridInterpolator(
        (data.variables["lat"][:], data.variables["lon"][:]), data.variables["orog"][:]
    )
    coords2 = coords.copy()
    coords2[:, 1] = np.where(coords2[:, 1] < 0, coords2[:, 1] + 360, coords2[:, 1])

    ele = interp(coords2)
    return ele
