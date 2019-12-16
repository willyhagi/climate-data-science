#!/usr/bin/env python


'''
 File Name: correl.py
 Description: Pearson Linear Correlation for Gridded Data and Climate Indices.
 Author: Willy Hagi
 E-mail: hagi.willy@gmail.com
 Python Version: 3.7.3
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

# --- pearson linear correlation (as simple as that)
pearson_r, p_values = corr(sst, nino34, dim='time', return_p=True)

# --- plotting
f, ax = plot.subplots(axwidth=6., tight=True,
                      proj='pcarree', proj_kw={'lon_0': 180},)
# format options
ax.format(land=False, coast=True, innerborders=True, borders=True,
          large='15px', labels=False,
          latlim=(31, -31), lonlim=(119, 291),
          geogridlinewidth=0,)

map1 = ax.contourf(dset['lon'], dset['lat'], pearson_r,
                   levels=np.arange(-1.0, 1.1, 0.1), cmap='Div', extend='both')
# still need to figure out how to plot p_values
# ax.contour(dset['lon'], dset['lat'], p_values,
#           colors=('k',),linestyles=('--','-'),
#           labels=True)

ax.colorbar(map1, loc='b', shrink=0.5, extendrect=True)

plt.show()
