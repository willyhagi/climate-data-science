#!/usr/bin/env python


'''
 File Name: correl.py
 Description: Pearson Linear Correlation for Gridded Data and Climate Indices.
 Observations: Statistical Significance is left for the user determination.
 Author: Willy Hagi
 E-mail: hagi.willy@gmail.com
 Python Version: 3.6
'''





import cartopy.crs          as ccrs
import cartopy.feature      as cf
import cartopy              as cartopy
import matplotlib.pyplot    as plt
import matplotlib.gridspec  as gridspec
import numpy                as np
import xarray               as xr
import scipy

from cartopy.mpl.ticker    import LongitudeFormatter, LatitudeFormatter
from scipy.signal          import detrend
from cartopy.util 		   import add_cyclic_point
from functions             import*





##----- read files
dset  =  xr.open_dataset('pacific_asst_1951_2015.nc')
sst   =  dset['sst'][:,:,:]
lat   =  dset['lat'][:]
lon   =  dset['lon'][:]


nino3 = np.loadtxt('nino3.asc')



##--- numpy array converting (useful for plotting)
lon  =  np.asarray(lon.lon.values)
lat  =  np.asarray(lat.lat.values)





##----------------------- CALCULATIONS

t = (sst['time'])  # number of months

##--- masking
asst_nan   =  Nan_calc(sst, t, lat, lon)
asst_res   =  asst_nan.reshaping()
val, asst_masked   =  asst_nan.masking_array()

grid_mask =  len(asst_masked[0,0:])
time      =  len(asst_masked[0:,0])



##--- detrending data
asst_masked  =  detrend(asst_masked, type='linear')
nino3        =  detrend(nino3, type='linear')



##--- pearson linear correlation
r, p = pearson(asst_masked, nino3, grid_mask)
r = r.reshape(1, grid_mask)


##--- t-student test
t1  =  t_test(r, 65, 2.306)     # number of years



##--- matrix reconstruction
naxis   =  [np.newaxis]

r_rec   =  rec_matrix(r, naxis, lat, lon, val)
t1_rec  =  rec_matrix(t1, naxis, lat, lon, val)





##----------------------- PLOTTING

plt.figure(figsize=(12,5))
inter = np.arange(-1.0, 1.2, 0.1)
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

plt.contour(lon, lat, t1_rec[0,:,:],interc,
            colors=('k',),linestyles=('--','-'),
            transform=ccrs.PlateCarree(),)


plt.contourf(lon, lat, r_rec[0,:,:], inter,
             transform=ccrs.PlateCarree(),
             cmap=plt.get_cmap('RdBu_r'))

plt.colorbar(ax=ax, shrink=0.5, orientation='horizontal')
plt.title(u'Correl [ASST, Niño 3]')

plt.tight_layout()
plt.show()





# A. M. D. G. #