#!/usr/bin/env python


'''
 File Name: annual_cycle.py
 Description: Gridded monthly average and anomaly values.
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

from calendar              import month_name
from functions             import* # from functions.py file
from cartopy.mpl.ticker    import LongitudeFormatter, LatitudeFormatter
from scipy.signal          import detrend





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

fclim  =  Demean(var, lat1, lat2, lon1, lon2, time1, time2)  # call Demean class
clim   =  fclim.average()                                    # monthly mean
anom   =  fclim.anomaly()                                    # anomalies
index  =  fclim.area_avg()                                   # index calculation



##----- DETRENDING (OPTIONAL)
''' Note: it is highly advisable to detrend you data at this point. It can save you from some unnecessary 
repetitions of lines of code later (when, for example, if you want to perform some EOF analysis), avoiding
any kind of conflict between Xarray and SciPy AND making your things run faster. The following lines 
are mostly from Dr. Nicholas Fauchereau's notebooks (http://nicolasfauchereau.github.io/climatecode/posts/).'''


##--- processing of composite matrix
anom = np.asarray(anom)
naxis  =  [np.newaxis]

# masking
asst_nan   =  Nan_calc(anom, time, lat, lon)
asst_res   =  asst_nan.reshaping()
val, anom  =  asst_nan.masking_array()

# detrend
anom = detrend(anom, type='linear')

# rec 
anom   =  rec_matrix(anom, time, lat, lon, val) 


d = {}
d['time'] = ('time',time)
d['lat'] = ('lat',lat)
d['lon'] = ('lon', lon)
d['sst'] = (['time','lat','lon'], anom)

dset_from_dict = xr.Dataset(d)
print (type(dset_from_dict))

dset_from_dict.to_netcdf('asstv5.nc')




#--- saving files
#anom = anom.to_dataset()
#anom.to_netcdf('aprp_teste.nc')            # save anomaly values to netcdf file
#np.savetxt('index.asc',index,fmt='%.10f')  # save index to text file





##----------------------- PLOTTING



interval  = np.arange(0,700,100)
proj      =  ccrs.PlateCarree()     # cartopy projection

fig, axes = plt.subplots(figsize=(10,12),nrows = 3, 
						 ncols = 4, subplot_kw={'projection': proj})

axes = axes.flatten()

for i, month in enumerate(range(1,13,1)):
    ax = axes[i]
#    ax.set_extent([280, 325, -50, 12], proj)                 # select lon/lat
    ax.add_feature(cf.BORDERS,zorder=2,edgecolor='k')        # borders
    ax.add_feature(cartopy.feature.COASTLINE)                # coastline
    x_lons = np.arange(-85,-15,20)
    y_lats = np.arange(-70,50,10)
    tick_fs = 16
    ax.set_xticks(x_lons, minor=False, crs=proj)
    ax.set_yticks(y_lats, minor=False, crs=proj)
    lon_formatter = LongitudeFormatter(zero_direction_label=True,
                                   number_format='.0f')
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.set_title(month_name[i+1])  # month titles
    data = ax.contourf(lon,lat,clim[i,:,:],interval,
                       transform = proj,
                       cmap = 'rainbow')
    fig.colorbar(data, ax=ax, shrink=0.5)

plt.tight_layout()
#plt.savefig('fig.eps',format='eps')  # save figure
plt.show()





# A. M. D. G. #
