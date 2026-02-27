import numpy as np
import matplotlib.pyplot as plt
from read_iasi_l1c import Read_Iasi_L1c  # Anu's code for reading iasi l1c data

import warnings
import matplotlib.collections
import matplotlib.ticker as mticker
import matplotlib as mpl
import math
from scipy.spatial import KDTree

"""Define some useful functions first."""


def decimal_degree_to_DMS(dec):
    negative = dec < 0
    dec = abs(dec)
    M, S = divmod(dec * 3600, 60)
    D, M = divmod(M, 60)
    if negative:
        if D < 0:
            D = -D
        elif M < 0:
            M = -M
        elif S < 0:
            S = -S
    return (int(D), int(M), int(S))


def lst_closest(lst, val):
    """Find closest to value in a list."""
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - val))]


def lst_closest_idx(lst, val):
    """Find index of closest to value in a list."""
    return min(range(len(lst)), key=lambda i: abs(lst[i] - val))


def get_channel_bbt(wno_spec, bbt_spec, wncenter=None, wnrange=None):
    """Get brightness temperature in a channel from a brightness temperature
    spectrum.
    wno_spec - wavenumbers list (x axis of spectrum)
    bbt_spec - brightness temperature spectrum (y axis of the spectrum)
    wncenter - tuple
               desired channel central wavenumber and width around that
               width units = central wn units
    wnrange - tuple, alternative way to specify the desired channel
              desired min and max wavenumber
              finds closest to each of these values
    """

    if wncenter == None and wnrange == None:
        return print("Please specify channel wavenumber or range.")
    elif wncenter != None and wnrange != None:
        return print("Please specift only channel wavenumber or range, not both.")
    elif wncenter != None and wnrange == None:
        if type(wncenter) == tuple:
            wmin = wncenter[0] - wncenter[1]
            wmax = wncenter[0] + wncenter[1]
            if wmin < wno_spec[0] or wmax > wno_spec[-1]:
                return print("Warning! Out of bounds! Your channel is too wide.")
        else:
            return print("Wncenter is not a tuple!")
    elif wncenter == None and wnrange != None:
        if type(wnrange) == tuple:
            wmin = wnrange[0]
            wmax = wnrange[1]
        else:
            return print("Wnrange is not a tuple!")
    else:
        return print("Something unexpected is wrong, check all inputs again.")

    wmin_closest_idx = lst_closest_idx(wno_spec, wmin)
    wmax_closest_idx = lst_closest_idx(wno_spec, wmax)

    if wmin_closest_idx == wmax_closest_idx:
        bbt = bbt_spec[wmin_closest_idx]
    else:
        if wmax_closest_idx + 1 > len(bbt_spec):
            wmax_closest_idx = -1
        bbt = sum(bbt_spec[wmin_closest_idx:wmax_closest_idx]) / len(
            bbt_spec[wmin_closest_idx:wmax_closest_idx]
        )
    return bbt


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

"""Define the desired channel/range for brightness temperature."""
wncenter = (833.3, 1)

"""Define some other things - symbols, constants, etc."""
deg_sym = "\u00b0"
min_sym = "'"
sec_sym = '"'

"""Read the IASI .nat file."""
l1c = Read_Iasi_L1c(
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

"""Calculate brightness temperature"""
bbts = []
for spec in l1c.spec:
    bbt = get_channel_bbt(l1c.wno, spec, wncenter=wncenter)
    bbts.append(bbt)

"""Plot spectrum."""
print(f"Attempting to plot brightness temperature histogram.")

plt.ion()
plt.figure(figsize=(11, 8.5))

ax = plt.subplot(1, 1, 1)
ax.hist(bbts, bins=100, edgecolor="k")

ax.set_xlabel("Brightness temperature(K)")
ax.set_ylabel("No. of occurences")
ax.set_title(f"Brightness temperature histogram")
# ax.legend()

# plt.savefig('test.pdf',dpi=1200)
# plt.tight_layout()
plt.show()

print("Succesful, there you go.")
