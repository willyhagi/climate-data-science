#!/usr/bin/env python


'''
 File Name: cartopy_plot.py
 Description: Introduction to Cartopy Plotting with Xarray.
 Observations: None.
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

from cartopy.mpl.ticker    import LongitudeFormatter, LatitudeFormatter





##----- read netcdf file
dset  =  xr.open_dataset('precip.mon.total.1x1.v7.nc')
var   =  dset['precip'][:,:,:]
lat   =  dset['lat'][:]
lon   =  dset['lon'][:]



##----- if you want to select a certain area (and time) over the globe
# south america
lat1  =  12 ; lat2  =  -50 ; lon1  =  270 ; lon2  =  330  


var  =  var.sel(lat=slice(lat1,lat2),
				lon=slice(lon1,lon2))



##----- numpy array converting (useful for plotting)
lon  =  var.sel(lon=slice(lon1,lon2))  
print (type(lon))
lon  =  np.asarray(lon.lon.values)
print (type(lon))

lat  =  var.sel(lat=slice(lat1,lat2)) 
lat  =  np.asarray(lat.lat.values)





##----------------------- PLOTTING
plt.figure(figsize=(8,4))
proj  =  ccrs.PlateCarree()
ax = plt.axes(projection=proj)
inter = np.arange(0, 800, 100)

ax.set_extent([250, 340, -50, 15], proj)             
ax.add_feature(cf.BORDERS)
ax.add_feature(cf.COASTLINE)
ax.coastlines(resolution='50m',color='black')
x_lons = np.arange(-105,-15,15)
y_lats = np.arange(-50,30,20)
tick_fs = 16
ax.set_xticks(x_lons, minor=False, crs=proj)
ax.set_yticks(y_lats, minor=False, crs=proj)
lon_formatter = LongitudeFormatter(zero_direction_label=True,
                                   number_format='.0f')
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

plt.contourf(lon, lat, var[0,:,:], inter,
             transform=proj,
             cmap=plt.get_cmap('RdBu'))
plt.colorbar(ax=ax, shrink=0.8)
plt.title(u'Precipitation over South America (mm/month)')

plt.show()





# A. M. D. G. #
