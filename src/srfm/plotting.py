"""Provides functions useful for plotting data that are not default in python.

- Name: plotting
- Parent package: srfm
- Author: Antonin Knizek
- Contributors:
- Date: 18 February 2025
"""


def map_extent(center_lon, center_lat, lon_span, lat_span):
    """Calculates cartopy map extent.

    Calculation performed from a known centrepoint coordinates + span on each side.
    Units decimal degrees.
    Remember that central_longitude of the projection is now your 0.

    Args:
        center_lon (float): Centrepoint longitude.
        center_lat (float): Centrepoint latitude.
        lon_span (float): Longitude span.
        lat_span (flaot): Latitude span.

    Returns:
        l (list of float): list containing [longitude minimum, longitude maximum,
            latitude minimum, latitude maximum].

    """
    lon_min = -lon_span
    lon_max = lon_span

    if center_lat - lat_span < -90:
        lat_min = center_lat - lat_span + 180
    else:
        lat_min = center_lat - lat_span

    if center_lat + lat_span > 90:
        loat_max = center_lat + lat_span - 180
    else:
        lat_max = center_lat + lat_span

    l = [lon_min, lon_max, lat_min, lat_max]
    return l
