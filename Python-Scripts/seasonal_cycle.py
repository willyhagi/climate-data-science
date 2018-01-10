#!/usr/bin/env python


'''
 File Name: seasonal_cycle.py
 Description: Seasonal analysis with Xarray.
 Observations: Your input data must be a three-dimensional netcdf file.
 Author: Willy Hagi
 E-mail: hagi.willy@gmail.com
 Python Version: 3.6
'''





import matplotlib.pyplot   as plt
import cartopy.crs         as ccrs
import cartopy.feature     as cf
import cartopy             as cartopy
import numpy               as np
import xarray              as xr
import datetime


from functions    		   import*
from cartopy.mpl.ticker    import LongitudeFormatter, LatitudeFormatter





##----- read netcdf file
dset  =  xr.open_dataset('precip.mon.total.1x1.v7.nc')
var   =  dset['precip'][:,:,:]
lat   =  dset['lat'][:]
lon   =  dset['lon'][:]



##----- if you want to select a certain area (and time) over the globe

# south america
lat1  =  12  
lat2  =  -50  
lon1  =  270  
lon2  =  330  


time1 = '1981-1-1' ; time2 = '2010-12-1'



##----- numpy array converting (useful for plotting)
lon  =  var.sel(lon=slice(lon1,lon2))  
lon  =  np.asarray(lon.lon.values)

lat  =  var.sel(lat=slice(lat1,lat2)) 
lat  =  np.asarray(lat.lat.values)





##----------------------- CALCULATIONS

#--- user-defined seasons
djf  =  Seasonal(var, 1, 2,  12 )   
djf  =  djf.season()

mam  =  Seasonal(var, 3, 4,  5  )  
mam  =  mam.season()

jja  =  Seasonal(var, 6, 7,  8  )  
jja  =  jja.season()

son  =  Seasonal(var, 9, 10, 11 )  
son  =  son.season()



#--- gridded seasonal averages
idjf  =  Demean(djf, lat1, lat2, lon1, lon2, time1, time2)
mdjf  =  idjf.average()     # you still have D-J-F here
mdjf  =  mdjf.mean(axis=0)  # now you have a seasonal average

imam  =  Demean(mam, lat1, lat2, lon1, lon2, time1, time2)
mmam  =  imam.average()
mmam  =  mmam.mean(axis=0)

ijja  =  Demean(jja, lat1, lat2, lon1, lon2, time1, time2)
mjja  =  ijja.average()
mjja  =  mjja.mean(axis=0)

ison  =  Demean(son, lat1, lat2, lon1, lon2, time1, time2)
mson  =  ison.average()
mson  =  mson.mean(axis=0)





##----------------------- PLOTTING
interval  =  np.arange(0,600,100)  # PRP range
proj      =  ccrs.PlateCarree()    # cartopy projection

x_lons = np.arange(-85,-15,20) 
y_lats = np.arange(-70,50,10)
tick_fs = 16
lon_formatter = LongitudeFormatter(zero_direction_label=True,
                                   number_format='.0f')
lat_formatter = LatitudeFormatter()
    

fig, axes = plt.subplots(figsize=(6,6),nrows = 2, ncols = 2,
                        subplot_kw={'projection': proj},
                        sharex=False)
ax = axes


# DJF
ax[0,0].add_feature(cf.BORDERS,zorder=2,edgecolor='k')
ax[0,0].add_feature(cartopy.feature.COASTLINE)
ax[0,0].set_xticks(x_lons, minor=False, crs=proj)
ax[0,0].set_yticks(y_lats, minor=False, crs=proj)
ax[0,0].xaxis.set_major_formatter(lon_formatter)
ax[0,0].yaxis.set_major_formatter(lat_formatter)

data1 = ax[0,0].contourf(lon, lat, mdjf,interval,
                         transform=proj,cmap='Blues')
fig.colorbar(data1, ax=ax[0,0], shrink=0.5)
ax[0,0].set_title('DJF')


# MAM
ax[0,1].add_feature(cf.BORDERS,zorder=2,edgecolor='k')
ax[0,1].add_feature(cartopy.feature.COASTLINE)
ax[0,1].set_xticks(x_lons, minor=False, crs=proj)
ax[0,1].set_yticks(y_lats, minor=False, crs=proj)
ax[0,1].xaxis.set_major_formatter(lon_formatter)
ax[0,1].yaxis.set_major_formatter(lat_formatter)
data2 = ax[0,1].contourf(lon, lat, mmam,interval,
                         transform=proj,cmap='Blues')
fig.colorbar(data1, ax=ax[0,1], shrink=0.5)
ax[0,1].set_title('MAM')


# JJA
ax[1,0].add_feature(cf.BORDERS,zorder=2,edgecolor='k')
ax[1,0].add_feature(cartopy.feature.COASTLINE)
ax[1,0].set_xticks(x_lons, minor=False, crs=proj)
ax[1,0].set_yticks(y_lats, minor=False, crs=proj)
ax[1,0].xaxis.set_major_formatter(lon_formatter)
ax[1,0].yaxis.set_major_formatter(lat_formatter)
data2 = ax[1,0].contourf(lon, lat, mjja,interval,
                         transform=proj,cmap='Blues')
fig.colorbar(data1, ax=ax[1,0], shrink=0.5)
ax[1,0].set_title('JJA')


# SON
ax[1,1].add_feature(cf.BORDERS,zorder=2,edgecolor='k')
ax[1,1].add_feature(cartopy.feature.COASTLINE)
ax[1,1].set_xticks(x_lons, minor=False, crs=proj)
ax[1,1].set_yticks(y_lats, minor=False, crs=proj)
ax[1,1].xaxis.set_major_formatter(lon_formatter)
ax[1,1].yaxis.set_major_formatter(lat_formatter)

data2 = ax[1,1].contourf(lon, lat, mson,interval,
                         transform=proj,cmap='Blues')
fig.colorbar(data1, ax=ax[1,1], shrink=0.5)
ax[1,1].set_title('SON')


plt.tight_layout()

plt.show()





# A. M. D. G #