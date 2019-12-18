#!/usr/bin/env python


'''
 File Name: hovmoller_plot.py.
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
# something around the amazon basin
dset = dset.sel(lat=slice(5, -20), lon=slice(280, 310),
                time=slice('1981-01-01', '2010-12-01'))

# --- time and lat/lon averaging
hovmoller = dset['precip'].groupby('time.month').mean('time')
hovmoller_lat = hovmoller.mean('lon')
hovmoller_lon = hovmoller.mean('lat')


# --- plotting
f, ax = plot.subplots(figsize=(7, 5), ncols=2, tight=True, sharex=False)
# format options
ax.format(small='15px', large='15px', abc=True, abcloc='ul')

months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
          'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']

ax[0].contourf(hovmoller_lat['lat'], hovmoller['month'], hovmoller_lat,
               levels=np.arange(0, 450, 50), cmap='Dusk', labels=True,
               extend='both')
ax[1].contourf(hovmoller_lon['lon'], hovmoller['month'], hovmoller_lon,
               levels=np.arange(0, 450, 50), cmap='Dusk', labels=True,
               extend='both')

ax[0].format(xformatter='deglat', ylocator='index',
             yformatter=months, xlabel='', ylabel='', title='')
ax[1].format(xformatter='deglon',
             xlabel='', ylabel='', title='')

f.save('hovmoller.jpeg', dpi=300)

plt.show()
