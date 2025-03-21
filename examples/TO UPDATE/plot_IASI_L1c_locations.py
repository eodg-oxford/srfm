import oxharp
import numpy as np
import matplotlib.pyplot as plt
from cartopy import crs as ccrs, feature as cfeature
from cartopy.mpl.ticker import (LatitudeFormatter, LongitudeFormatter,
                                LatitudeLocator,LongitudeLocator)
import warnings
import matplotlib.collections
import matplotlib.ticker as mticker

earth_circumference = 40075 #Earth circumference, units km
iasi_pix = 12 #iasi pixel diamater on ground at nadir, units km
HTHH_lat = -20.545 #units decimal degrees
HTHH_lon = -175.3925 #units decimal degrees

###READ IASI NAT FILE###
fl = '/network/aopp/apres/IASI-C/2022/01/IASI_xxx_1C_M03_20220119223255Z_20220120001455Z_N_O_20220120001322Z.nat'
    #iasi filename, contains one whole orbit of data
avhrr = False
    #bool, True = include AVHRR cluster analysis data.
bright = True
    #bool, True for brightness temp spectra, False for radiance spectra
chkqal = (True,True,True)
    #bool, True = Exclude data with bad quality flag for each band, 
    #False = include, (bands are 645-1210, 1210-2000, 2000-2760 cm-1). 
cldlim = None
    #tuple, min,max cloud cover percentages for inclusion, (0,0)=only cloud-free
latlim = None
    #tuple, range  min,max latitudes (range -90:90) for inclusion
lndlim = None
    #tuple, min,max land (cf.ocean) percentages for inclusion, (0,0)=only ocean
loc_only = True
    #bool, Set True if only location data (without spectra) are required
lonlim = None
    #tuple, min,max longitudes (range -180:180) for inclusion, default=None
    #lonmin,lonmax are actually east,west boundaries so, for example, 
    #(-170,170) would include everything apart from the 20deg lon band 
    #around +/-180 deg, whereas (170,-170) would only include this 
    #20deg band.
    #The same actually applies to all other limits, but is probably 
    #only useful for longitude
mph_only = False
    #bool, Set True if just the Main Product Header data is required 
    #(as self.mph)
szalim = None
    #tuple, min,max solar zenith angle (0,90=daytime, 90-180=nighttime).
wnolim = (700,1400)
    #tuple, min,max wavenumber [cm-1] limits for extracted spectra.
zenlim = None
    #tuple, min,max satellite zenith angle (0=sat overhead). Default=None.

l1c = oxharp.readers.read_iasi_l1c.Read_Iasi_L1c(fl,
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
                    zenlim=zenlim
                   )
if l1c.fail:
    print(l1c.errmsg)
    stop

###PLOT DATA LOCATIONS###
plt.ion()
plt.figure(figsize=(11,8.5))
projPC = ccrs.PlateCarree(central_longitude=HTHH_lon)
ax = plt.subplot(1,1,1,projection=projPC)

ax.coastlines()
#ax.add_feature(cfeature.OCEAN)
#do not use ocean colour (facecolor) for readabilty, do not use borders for political reasons
meas_pt_size=iasi_pix*360/earth_circumference #sets point size to represent area iasi sees on ground at nadir 
ax.scatter(l1c.lon,
           l1c.lat,
           transform=ccrs.PlateCarree(),
           c='darkorange',
           s=meas_pt_size) #here the transform argument says which projection are our data in, also PlateCarree equals Cartesian, does not work on zoom
#meas_pts = [plt.Circle((loni,lati),radius=meas_pt_size) for loni,lati in zip(l1c.lon,l1c.lat)]
#p = matplotlib.collections.PathCollection(meas_pts)
#ax.add_collection(p)
ax.scatter(HTHH_lon,
           HTHH_lat,
           c='r',
           s=200,
           transform=ccrs.PlateCarree())

#set map extent
extent=oxharp.plotting.map_extent(center_lon=HTHH_lon,
                         center_lat=HTHH_lat,
                         lon_span=50,
                         lat_span=50
                         )

ax.set_extent(extent,
              crs=ccrs.PlateCarree(central_longitude=HTHH_lon)
              )
              
'''    
Extent cannot span longitude more than 360.
If you specify central longitude in Plate Carree as 0 or do not specify (0 is inferred), you get the usual lat/lon coordinate system
If you specify a different central longitude, then this becomes your new zero and your span is still -180,180, but relative to this new zero.
So if you set central_location=18.068 in the PlateCarree coordidante system, your map will be centered around Stockholm and London will now be at -18.068.
In this way you can make plots that cross the dateline.
Make sure you set the central location when setting up the GeoAxes and use the same when specyfing the extent
'''

#set ticks manually:
#ax.set_xticks([-180,-177.5,-175,-172.5,-170,-167.5,177.5],crs=ccrs.PlateCarree())
#lon_formatter = LongitudeFormatter(zero_direction_label=True)
#ax.xaxis.set_major_formatter(lon_formatter)
#ax.set_yticks([HTHH_lat-8,HTHH_lat-4,HTHH_lat,HTHH_lat+4,HTHH_lat+8],crs=ccrs.PlateCarree())
#lat_formatter = LatitudeFormatter()
#ax.yaxis.set_major_formatter(lat_formatter)
#'''Overshooting plot boundaries changes them
#ticks must be in increasing order and within -180,180 range,
#so if crossing the dateline, reorder the list accortingly.
#'''

#set ticks automatically:
gl = ax.gridlines(draw_labels=['bottom','left'], #draws labels
                  dms=False,
                 )

#sets tick spacing:
gl.xlocator = LongitudeLocator(nbins=180) #for some reason sets the spacing as 360/nbins, if not specified, creates ticks every 36 deg, which is useless for zoomed plots
gl.ylocator = LatitudeLocator()
gl.xlines=False
gl.ylines=False

#plt.savefig('test.pdf',dpi=1200)
#plt.tight_layout()
plt.show()

