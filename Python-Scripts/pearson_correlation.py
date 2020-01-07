#!/usr/bin/env python


'''
 File Name: correl.py
 Description: Pearson Linear Correlation for Gridded Data and Climate Indices.
 Author: Willy Hagi
 E-mail: hagi.willy@gmail.com
 Python Version: 3.6.7
'''


import numpy as np
import xarray as xr
import proplot as plot
import matplotlib.pyplot as plt

from esmtools.stats import*


# --- read netcdf file
dset = xr.open_dataset('asstdt_pacific.nc')

# --- select djf months
sst = dset['sst'].sel(time=np.in1d(dset['time.month'], [1, 2, 12]))
print(sst)

# --- make niño 3.4 index
nino34 = sst.sel(lat=slice(5, -5), lon=slice(360 - 170, 360 - 120))
nino34 = nino34.mean(dim=('lat', 'lon'))

# --- pearson linear correlation
pearson_r, p_values = corr(sst, nino34, dim='time', return_p=True)

# --- plotting
fig, ax = plot.subplots(axwidth=6., tight=True,
                      proj='pcarree', proj_kw={'lon_0': 180},)
# format options
ax.format(land=False, coast=True, innerborders=True, borders=True,
          large='15px', labels=True,
          latlim=(31, -31), lonlim=(119, 291),
          lonlines=plot.arange(130, 280, 20),
          geogridlinewidth=0,)

# plot correlation values
map1 = ax.contourf(dset['lon'], dset['lat'], pearson_r,
                   levels=50, cmap='ColdHot', extend='both')
# plot p_values
ax.contourf(dset['lon'], dset['lat'], p_values,
            levels=np.arange(0, 0.05, 0.01), hatches=['....'], alpha=0)
# colorbar
ax.colorbar(map1, loc='b', shrink=0.5, extendrect=True)

ax.format(title='Correlation between Niño 3.4 Index and ASST')

plt.show()

