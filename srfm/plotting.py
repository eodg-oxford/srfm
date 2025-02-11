def map_extent(center_lon, center_lat, lon_span, lat_span):
    """calculates cartopy map extent from a known centrepoint
    coordinates + span on each side.
    Units decimal degrees.
    Remember that central_longitude of the projection is now your 0.
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
    return [lon_min, lon_max, lat_min, lat_max]
