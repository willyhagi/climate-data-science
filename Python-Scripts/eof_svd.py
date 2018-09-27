#!/usr/bin/env python


'''
 File Name: eof_svd.py
 Description: Empirical Orthogonal Functions (EOF) with Singular Value Decomposition (SVD).
 Observations: Your input data must be a three-dimensional netcdf file
 / It's much better if your input data is already detrended, normalized and so on...
 / Statistical Significance is left for user determination.
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
import pandas               as pd
import xarray               as xr
import scipy
import string

from cartopy.mpl.ticker     import LongitudeFormatter, LatitudeFormatter
from calendar               import month_name
from scipy.signal           import detrend
from functions              import*
from scipy.stats            import pearsonr





##----- read netcdf file
dset  =  xr.open_dataset('pacific_asst_1951_2015.nc')
sst   =  dset['sst'][:,:,:]
lat   =  dset['lat'][:]
lon   =  dset['lon'][:]





##----------------------- CALCULATIONS


##--- month selecting (optional)
djf  =  sst.sel(time=np.in1d(sst['time.month'], [1,2,12]))



##--- some useful numpy arrays 
tseason  =  np.asarray(djf['time'])  # months array
p_lat    =  np.asarray(djf['lat'])   # latitude array
p_lon    =  np.asarray(djf['lon'])   # longitude array
time     =  len(tseason)             # time array 



##--- pre-processing of anomaly matrix

# masking to use detrend function
asst_nan   =  Nan_calc(djf, tseason, lat, lon)
asst_res   =  asst_nan.reshaping()
val, asst_masked   =  asst_nan.masking_array()

# linear detrending
asstdt   =  detrend(asst_masked, type='linear') 

# normalization of s-mode input matrix
F_matrix  =  asstdt / np.std(asstdt, axis=0)     


time       =  len(F_matrix[0:,0])
grid_mask  =  len(F_matrix[0,0:])



##--- EOF standard calculation and normalization 
U, V, g, lmbd, PC, EOF  =  eof_svd(F_matrix)

PCn, EOFn  =  eof_norm(PC, V.transpose(), lmbd)



##----- pearson correlation of principal components and ASST
PC1   =  PCn[0,:] # first PC
PC2   =  PCn[1,:] # second PC

#np.savetxt('PC1_V5.asc',PC1,fmt='%.10f')
#np.savetxt('PC2_V5.asc', PC2, fmt='%.10f')

r, p    =  pearson(asst_masked, PC1, grid_mask)
r       =  r.reshape(1, grid_mask)

r2, p2  =  pearson(asst_masked, PC2, grid_mask)
r2      =  r2.reshape(1, grid_mask)



##--- t-student test (mostly left for user determination)
t1  =  t_test(r, 47, 2.021)
t2  =  t_test(r2, 47, 2.021)


#--- matrix reconstructions
EOF_rec  =  rec_matrix(EOFn.transpose(), tseason, p_lat, p_lon, val)


naxis   =  [np.newaxis]
r_rec   =  rec_matrix(r, naxis, p_lat, p_lon, val)
r2_rec  =  rec_matrix(r2, naxis, p_lat, p_lon, val)

t1_rec  =  rec_matrix(t1, naxis, lat, lon, val)
t2_rec  =  rec_matrix(t2, naxis, lat, lon, val)



#----- amount of explained variance
var_exp = []
for i in range(1,len(np.arange(1,time+2,1)),1):
    var_exp.append((lmbd[i-1] / np.sum(lmbd)) * 100.)

var_exp = np.around(var_exp,decimals=2)

#np.savetxt('var_exp.asc', var_exp, fmt='%.10f')
#np.savetxt('lmbd.asc', lmbd, fmt='%.10f')





##----------------------- PLOTTING
n      =  len(tseason)
dt     =  (1/3)
time   =  (np.arange(n))*dt + 1951.0  

var_t  = np.arange(1, len(tseason)+1, 1)



##--- first and second principal components
fig, axes = plt.subplots(figsize=(14,9),nrows = 2, ncols = 1,)
ax = axes.flatten()

ax[0].plot(time, PC1 / np.std(PC1), 'k')
ax[0].fill_between(time, PC1/np.std(PC1), 0, where=PC1/np.std(PC1)>0, facecolor='red')
ax[0].fill_between(time, PC1/np.std(PC1), 0, where=PC1/np.std(PC1)<0, facecolor='blue')
ax[0].set_xticks(range(1950, 2020, 5))
ax[0].set_xlabel('Years', fontsize=14)
ax[0].set_ylabel(u'Amplitude', fontsize=14)
ax[0].set_title(u'a) CP2',
				   fontsize=16)

ax[1].plot(time, PC2 / np.std(PC2), 'k')
ax[1].fill_between(time, PC2/np.std(PC2), 0, where=PC2/np.std(PC2)>0, facecolor='red')
ax[1].fill_between(time, PC2/np.std(PC2), 0, where=PC2/np.std(PC2)<0, facecolor='blue')
ax[1].set_xticks(range(1950, 2020, 5))
#ax[1].set_yticks(range(-3,4, 1))
ax[1].set_xlabel('Years', fontsize=14)
ax[1].set_ylabel(u'Amplitude', fontsize=14)
ax[1].set_title(u'b) PC2',
				   fontsize=16)


plt.tight_layout(pad=4.)


plt.show()



#--- eigenvalue spectrum
plt.plot(var_t[0:10], var_exp[0:10])
plt.plot(var_t[0:10], var_exp[0:10], 'ro')
plt.title('Eigenvalue Spectrum')
plt.ylabel('Amount of Explained Variance')
plt.xlabel('EOF modes')
plt.xticks(range(1, 11, 1))

plt.tight_layout()
plt.show()



#---- correlation fields
inter   =  np.arange(-1.0, 1.1, 0.1)
interc  =  [-1.0, 1.0]
proj    =  ccrs.PlateCarree(central_longitude=180.)
proj2   =  ccrs.PlateCarree()

fig, ax = plt.subplots(figsize=(8,6),nrows=2, ncols=1,subplot_kw={'projection': proj})
ax = ax.flatten()

# a)
ax[0].add_feature(cf.BORDERS)
ax[0].add_feature(cf.COASTLINE)
ax[0].coastlines(resolution='50m',color='black')
x_lons  = np.arange(-180,360,60)
y_lats  = np.arange(-40,40,20)
tick_fs = 16
lon_formatter = LongitudeFormatter(zero_direction_label=True,
                                   number_format='.0f')
lat_formatter = LatitudeFormatter()
ax[0].xaxis.set_major_formatter(lon_formatter)
ax[0].yaxis.set_major_formatter(lat_formatter)
ax[0].set_xticks(x_lons, minor=False, crs=proj)
ax[0].set_yticks(y_lats, minor=False, crs=proj)
ax[0].set_title(u'a) CORREL [PC1, ASST]')

ts1 = ax[0].contour(lon, lat, t1_rec[0,:,:],interc,
                         colors=('k',),linestyles=('--','-'),
                         transform=ccrs.PlateCarree(),)

field1 = ax[0].contourf(lon, lat, r_rec[0,:,:], inter,
             transform=proj2,
             cmap='RdBu_r')  # RdBu_r / RdBu

fig.colorbar(field1,ax=ax[0], shrink=0.3, orientation='horizontal')

# b)
ax[1].add_feature(cf.BORDERS)
ax[1].add_feature(cf.COASTLINE)
ax[1].coastlines(resolution='50m',color='black')
y_lats  = np.arange(-40,40,10)
x_lons  = np.arange(-70,360,20)
tick_fs = 16
lon_formatter = LongitudeFormatter(zero_direction_label=True,
                                   number_format='.0f')
lat_formatter = LatitudeFormatter()
ax[1].xaxis.set_major_formatter(lon_formatter)
ax[1].yaxis.set_major_formatter(lat_formatter)
ax[1].set_xticks(x_lons, minor=False, crs=proj)
ax[1].set_yticks(y_lats, minor=False, crs=proj)
ax[1].set_title(u'b) CORREL [PC2, ASST]')

ts2 = ax[1].contour(lon, lat, t2_rec[0,:,:],interc,
                         colors=('k',),linestyles=('--','-'),
                         transform=ccrs.PlateCarree(),)

field2 = ax[1].contourf(lon, lat, r2_rec[0,:,:],inter ,
             transform=proj2,
             cmap='RdBu_r')  # RdBu_r / RdBu

fig.colorbar(field2, ax=ax[1], shrink=0.3, orientation='horizontal')  # barra de cores


plt.tight_layout()
plt.show()



#--- first 2 EOF modes
proj      =  ccrs.PlateCarree(central_longitude=180.)
proj2     =  ccrs.PlateCarree()


fig, axes = plt.subplots(figsize=(8,6),nrows = 2, ncols = 1,
                         subplot_kw={'projection': proj})
axes = axes.flatten()


for i, month in enumerate(range(1,3,1)):
    ax = axes[i]
    ax.add_feature(cf.BORDERS,zorder=2,edgecolor='k')           
    ax.add_feature(cartopy.feature.COASTLINE)                   
#    ax.add_feature(cf.LAND,color='grey')                      
    y_lats  = np.arange(-40,40,10)
    x_lons  = np.arange(-70,360,20)
    tick_fs = 16
    lon_formatter = LongitudeFormatter(zero_direction_label=True,
                                       number_format='.0f')
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.set_xticks(x_lons, minor=False, crs=proj)
    ax.set_yticks(y_lats, minor=False, crs=proj)
    ax.set_title('%s) EOF%s - %s%%' %(string.ascii_lowercase[i], i+1,var_exp[i])) # var_exp[i]
    data = ax.contourf(lon,lat,EOF_rec[i,:,:],
                       60,
                       transform = proj2,
                       cmap = 'RdBu_r')
    fig.colorbar(data, ax=ax, shrink=0.65)

plt.tight_layout()
plt.show()





# A. M. D. G. #
