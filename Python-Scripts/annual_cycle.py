#!/usr/bin/env python


'''
 File Name: annual_cycle.py
 Description: Gridded monthly average and anomaly values.
 Observations: Your input data must be a three-dimensional netcdf file.
 Author: Willy Hagi
 E-mail: hagi.willy@gmail.com
 Python Version: 3.6
'''


import numpy as np
import xarray as xr
import proplot as plot
import matplotlib.pyplot as plt

from esmtools.stats import*
from calendar import month_name


##--- read netcdf file
dset = xr.open_dataset('precip.mon.total.v2018.nc')
#--- select an area and time (optional)
dset = dset.sel(lat=slice(15, -60), lon=slice(270, 330))

##--- monthly climatology and anomalies with respect to the 1981-2010 period
base_period = dset['precip'].sel(time=slice('1981-1-1', '2010-12-1'))
wrt = base_period.groupby('time.month').mean('time')
anomalies = dset['precip'].groupby('time.month') - wrt


##--- detrend and normalize by the standard deviation
anom_dt = rm_poly(anomalies, 1, dim='time')
anom_norm = anom_dt / base_period.std('time')

##--- save files (optional)
#anom_dt.to_netcdf('asstdt.nc')
#anom_norm.to_netcdf('asstdt_norm.nc')


##--- plotting
f, ax = plot.subplots(axwidth=1.5, nrows=3, ncols=4, tight=True,
                      proj='pcarree', proj_kw={'lon_0': 0},)

ax.format(land=False, coast=True, innerborders=True, borders=True,
          large='15px', labels=True,
          latlim=(16, -51), lonlim=(269, 329),
          geogridlinewidth=0,
          abcloc='ur',
          )

# don't know why the month title is not working
for i, month in enumerate(range(1,13,1)):
    data = ax.contourf(wrt['lon'], wrt['lat'], wrt[i,:,:],
                       levels = np.arange(0, 550, 50),
                       cmap = 'BuPu',
                       extend='both')
    ax.format(title='{:s}'.format(month_name[i+1]))

plt.show()
