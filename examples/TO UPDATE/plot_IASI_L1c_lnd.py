import oxharp
import numpy as np
import matplotlib.pyplot as plt
from cartopy import crs as ccrs, feature as cfeature
from cartopy.mpl.ticker import (LatitudeFormatter, LongitudeFormatter,
                                LatitudeLocator,LongitudeLocator)
import matplotlib.collections
import matplotlib.ticker as mticker
import matplotlib as mpl

earth_circumference = 40075 #Earth circumference, units km
iasi_pix = 12 #iasi pixel diamater on ground at nadir, units km
HTHH_lat = -20.545
HTHH_lon = -175.3925

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
cldlim = (0,100)
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

###PLOT CLOUD COVER AS COLORMAP###
plt.ion()
plt.figure(figsize=(11,8.5))
projPC = ccrs.PlateCarree(central_longitude=HTHH_lon)
ax = plt.subplot(1,1,1,projection=projPC)

ax.coastlines()

meas_pt_size=iasi_pix*360/earth_circumference #sets point size to represent area iasi sees on ground at nadir 

##plor cloud cover
cmap='viridis'
ax.scatter(l1c.lon,
           l1c.lat,
           c=l1c.lnd,
           transform=ccrs.PlateCarree(),
           cmap=cmap,
           alpha=0.5
          )

#manually set colorbar (smh the object above is not mappable??)
norm = mpl.colors.Normalize(vmin=l1c.cld.min(), vmax=l1c.cld.max())
sm = plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
plt.colorbar(sm,ax=ax,label='Land cover (%)')

ax.scatter(HTHH_lon,
           HTHH_lat,
           c='r',
           s=200,
           transform=ccrs.PlateCarree()) #https://www.fcc.gov/media/radio/dms-decimal
       
extent=oxharp.plotting.map_extent(center_lon=HTHH_lon,
                         center_lat=HTHH_lat,
                         lon_span=50,
                         lat_span=50
                         )
ax.set_extent(extent,
              crs=ccrs.PlateCarree(central_longitude=HTHH_lon)
              )

#set ticks automatically:
gl = ax.gridlines(draw_labels=['bottom','left'], #draws labels
                  dms=False,
                 )

#sets tick spacing:
#gl.xlocator = mticker.MultipleLocator(1) 
gl.xlocator = LongitudeLocator(nbins=180) #for some reason sets the spacing as 360/nbins, if not specified, creates ticks every 36 deg, which is useless for zoomed plots
gl.ylocator = LatitudeLocator()
gl.xlines=False
gl.ylines=False

ax.set_xlabel('Longitude')
ax.set_ylabel('Lattitude')
ax.set_title('Percentage land cover')

#plt.savefig('test.pdf',dpi=1200)
#plt.tight_layout()
plt.show()

