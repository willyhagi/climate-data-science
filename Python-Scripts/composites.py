#!/usr/bin/env python


'''
 File Name: composites.py
 Description: El Niño composites.
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


''' the selection of the El Niño years are taken from here:
	http://ggweather.com/enso/oni.htm '''

# onset D(0) strong el niño years
onset_elnino = [1957, 1965, 1972, 1987, 1991]
# end JF(+1) strong el niño years
end_elnino = [1958, 1966, 1973, 1988, 1992]


# --- D(0)JF(+1) composites
# all the december months
dcb = dset['sst'].sel(time=np.in1d(dset['sst']['time.month'], [12]))
# D(0)
dcb0 = dcb.sel(time=np.in1d(dcb['time.year'], onset_elnino))
# JF(+1)
jf = dset['sst'].sel(time=np.in1d(dset['sst']['time.month'], [1, 2]))
jf1 = jf.sel(time=np.in1d(jf['time.year'], end_elnino))
# merge all months together
djf = xr.merge([dcb0, jf1])
# seasonal mean
djf = djf.groupby('time.season').mean('time')
print(djf)


# --- plotting
f, ax = plot.subplots(axwidth=5., tight=True,
                      proj='pcarree', proj_kw={'lon_0': 180},)
# format options
ax.format(land=False, coast=True, innerborders=True, borders=True,
          large='15px', labels=False,
          latlim=(31, -31), lonlim=(119, 291),
          geogridlinewidth=0)

map1 = ax.contourf(dset['lon'], dset['lat'], djf['sst'][0, :, :],
                   levels=np.arange(-1.5, 1.6, 0.1), cmap='Div', extend='both')

ax.colorbar(map1, loc='b', shrink=0.5, extendrect=True)

plt.show()
