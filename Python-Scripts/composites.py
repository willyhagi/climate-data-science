#!/usr/bin/env python


'''
 File Name: composites.py
 Description: El Niño composites.
 Observations: Statistical Significance is left for the user determination.
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

from calendar              import month_name
from functions             import*
from cartopy.util 		   import add_cyclic_point
from cartopy.mpl.ticker    import LongitudeFormatter, LatitudeFormatter
from scipy.signal          import detrend





##----- read netcdf file
dset  =  xr.open_dataset('pacific_asst_1951_2015.nc')
var   =  dset['sst'][:,:,:]
lat   =  dset['lat'][:]
lon   =  dset['lon'][:]



##--- numpy array converting (useful for plotting)
lon  =  np.asarray(lon.lon.values)
lat  =  np.asarray(lat.lat.values)





##----------------------- CALCULATIONS

''' the composite method comes from the following URL:
 	https://github.com/pydata/xarray/issues/557 '''


#--- select years
''' the selection of the ENSO years are taken from here:
	http://ggweather.com/enso/oni.htm '''


#--- composite years

# onset D(0) strong el niño years
onset_stenyr = [1957, 1965, 1972, 1987, 1991]
# end JF(+1) strong el niño years
end_stenyr   = [1958, 1966, 1973, 1988, 1992]




#--- D(0)JF(+1) composites
# all the decembers
dcb  =  var.sel(time=np.in1d(var['time.month'], [12])) 
# D(0)
dcb0 =  dcb.sel(time=np.in1d(dcb['time.year'], onset_stenyr))
# JF(+1)
jf   =  var.sel(time=np.in1d(var['time.month'], [1,2]))
jf1  =  jf.sel(time=np.in1d(jf['time.year'], end_stenyr))  


en_djf  = dcb0.mean('time') + jf1.mean('time') 
print (en_djf.shape)


##--- processing of composite matrix
naxis  =  [np.newaxis]

# masking
asst_nan   =  Nan_calc(en_djf, naxis, lat, lon)
asst_res   =  asst_nan.reshaping()
val, en_djf  =  asst_nan.masking_array()

grid_mask =  len(en_djf[0,0:])
time      =  len(en_djf[0:,0])



##--- detrending data
en_djf  =  detrend(en_djf, type='linear')



##--- t-student test
t_endjf = t_comp(en_djf, 5, 2.571)



##--- matrix reconstructions
ren_djf  =  rec_matrix(en_djf, naxis, lat, lon, val)
rt_endjf  =  rec_matrix(t_endjf, naxis, lat, lon, val)





##----------------------- PLOTTING
plt.figure(figsize=(12,5))
inter = np.arange(-4.0, 4.5, 0.5)
interc = [-1.0, 1.0]

proj  =  ccrs.PlateCarree(central_longitude=180.)
ax = plt.axes(projection=proj)
y_lats  = np.arange(-40,40,10)
x_lons  = np.arange(-70,360,20)
lon_formatter = LongitudeFormatter(zero_direction_label=True,
                                   number_format='.0f')
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
tick_fs = 16
ax.set_xticks(x_lons, minor=False, crs=proj)
ax.set_yticks(y_lats, minor=False, crs=proj)
ax.add_feature(cf.LAND,color='grey')
ax.add_feature(cf.BORDERS)
ax.add_feature(cf.COASTLINE)
ax.coastlines(resolution='50m',color='black')

plt.contour(lon, lat, rt_endjf[0,:,:], interc,
            colors=('k',),linestyles=('--','-'),
            transform=ccrs.PlateCarree(),)


plt.contourf(lon, lat, ren_djf[0,:,:], inter,
             transform=ccrs.PlateCarree(),
             cmap=plt.get_cmap('RdBu_r'))

plt.colorbar(ax=ax, shrink=0.5, orientation='horizontal')
plt.title(u'DJF Strong El Niño Composites')

plt.tight_layout()
plt.show()





# A. M. D. G. #