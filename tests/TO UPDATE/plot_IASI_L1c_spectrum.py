# import oxharp
import numpy as np
import matplotlib.pyplot as plt
import warnings
import matplotlib.collections
import matplotlib.ticker as mticker
import matplotlib as mpl
import math
from scipy.spatial import KDTree
from srfm import *

"""Define IASI .nat file name and parameters to be read from it."""
fl = "/network/aopp/apres/IASI-C/2022/01/IASI_xxx_1C_M03_20220119223255Z_20220120001455Z_N_O_20220120001322Z.nat"

avhrr = False
# bool, True = include AVHRR cluster analysis data.
bright = True
# bool, True for brightness temp spectra, False for radiance spectra
chkqal = (True, True, True)
# bool, True = Exclude data with bad quality flag for each band,
# False = include, (bands are 645-1210, 1210-2000, 2000-2760 cm-1).
cldlim = None
# tuple, min,max cloud cover percentages for inclusion, (0,0)=only cloud-free
latlim = None
# tuple, range  min,max latitudes (range -90:90) for inclusion
lndlim = None
# tuple, min,max land (cf.ocean) percentages for inclusion, (0,0)=only ocean
loc_only = False
# bool, Set True if only location data (without spectra) are required
lonlim = None
# tuple, min,max longitudes (range -180:180) for inclusion, default=None
# lonmin,lonmax are actually east,west boundaries so, for example,
# (-170,170) would include everything apart from the 20deg lon band
# around +/-180 deg, whereas (170,-170) would only include this
# 20deg band.
# The same actually applies to all other limits, but is probably
# only useful for longitude
mph_only = False
# bool, Set True if just the Main Product Header data is required
# (as self.mph)
szalim = None
# tuple, min,max solar zenith angle (0,90=daytime, 90-180=nighttime).
wnolim = (700, 1400)
# tuple, min,max wavenumber [cm-1] limits for extracted spectra.
zenlim = None
# tuple, min,max satellite zenith angle (0=sat overhead). Default=None.

"""Define longitude and lattitude you want to plot at.
   The code will find the nearest measurement to your specified location."""
desired_lon = 70
desired_lat = 70

"""Define some other things - symbols, constants, etc."""
deg_sym = "\u00b0"
min_sym = "'"
sec_sym = '"'

"""Read the IASI .nat file."""
l1c = readers.read_iasi_l1c.Read_Iasi_L1c(
    fl,
    avhrr=avhrr,
    bright=bright,
    chkqal=chkqal,
    cldlim=cldlim,
    latlim=latlim,
    lndlim=lndlim,
    loc_only=loc_only,
    lonlim=lonlim,
    mph_only=mph_only,
    szalim=szalim,
    wnolim=wnolim,
    zenlim=zenlim,
)
if l1c.fail:
    print(l1c.errmsg)
    stop

"""Find the closest measurement to your specified location."""
iloc = utilities.closest(l1c.lon, desired_lon, l1c.lat, desired_lat)

lon = float(l1c.lon[iloc])
lat = float(l1c.lat[iloc])

"""Convert decimal degree to DMS and prepare strings to print out."""
lon_DMS = units.decimal_degree_to_DMS(lon)
lat_DMS = units.decimal_degree_to_DMS(lat)

print_lon_DMS = (
    f"{lon_DMS[0]}{deg_sym}{lon_DMS[1]}{min_sym}{lon_DMS[2]}{sec_sym}"
    + str(["W" if lon_DMS[0] < 0 else "E"][0])
)
print_lat_DMS = (
    f"{lat_DMS[0]}{deg_sym}{lat_DMS[1]}{min_sym}{lat_DMS[2]}{sec_sym}"
    + str(["N" if lat_DMS[0] > 0 else "S"][0])
)

print("Nearest measurement point found at {print_lon_DMS}, {print_lat_DMS}.")

"""Plot spectrum."""
print(
    f"Attempting to plot brightness temperature spectrum at {print_lon_DMS}, {print_lat_DMS}"
)

plt.ion()
plt.figure(figsize=(11, 8.5))

ax = plt.subplot(1, 1, 1)
ax.plot(
    l1c.get_spec(iloc, wnolim=wnolim, bright=bright)[1],
    l1c.get_spec(iloc, wnolim=wnolim, bright=bright)[0],
    label=f"{print_lon_DMS}, {print_lat_DMS}",
)

ax.set_xlabel(r"Wavenumber (cm$^{-1}$)")
ax.set_ylabel("Brightness temperature (K)")
ax.set_title(f"Brightness temperature spectrum")
ax.legend()

# plt.savefig('test.pdf',dpi=1200)
# plt.tight_layout()
plt.show()
