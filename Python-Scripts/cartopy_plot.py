#!/usr/bin/env python


'''
 File Name: cartopy_plot.py
 Description: Introduction to Cartopy Plotting with Xarray.
 Observations: Now with Proplot!
 Author: Willy Hagi
 E-mail: hagi.willy@gmail.com
 Python Version: 3.7.3
'''


import numpy as np
import xarray as xr
import proplot as plot
import matplotlib.pyplot as plt


# --- read netcdf file
dset = xr.open_dataset('precip.mon.total.v2018.nc')
# --- select an area and time (optional)
dset = dset.sel(lat=slice(15, -60), lon=slice(270, 330))


#--- plotting
f, ax = plot.subplots(axwidth=4.5, tight=True,
                      proj='pcarree', proj_kw={'lon_0': 0},)
# format options
ax.format(land=False, coast=True, innerborders=True, borders=True,
          large='15px', labels=True,
          latlim=(16, -51), lonlim=(269, 329),
          geogridlinewidth=0,
          abcloc='ur',)

map1 = ax.contourf(dset['lon'], dset['lat'], dset['precip'][0, :, :],
                   levels=np.arange(0, 550, 50), cmap='BuPu', extend='both')

ax.colorbar(map1, loc='b', shrink=0.5, extendrect=True)

plt.show()
