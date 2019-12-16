#!/usr/bin/env python


'''
 File Name: eof.py
 Description: Empirical Orthogonal Functions (EOF).
 Author: Willy Hagi
 E-mail: hagi.willy@gmail.com
 Python Version: 3.7.3
'''


import numpy as np
import xarray as xr
import proplot as plot
import matplotlib.pyplot as plt

from esmtools.stats import*
from eofs.xarray import Eof


# --- read netcdf file
dset = xr.open_dataset('asstdt_pacific.nc')

# --- select djf months
sst = dset['sst'].sel(time=np.in1d(dset['time.month'], [1, 2, 12]))

# --- square-root of cosine of latitude weights
coslat = np.cos(np.deg2rad(sst.coords['lat'].values))
wgts = np.sqrt(coslat)[..., np.newaxis]
# --- eof solver
solver = Eof(sst, weights=wgts)
# --- eof results
eofs = solver.eofsAsCorrelation(neofs=2)
pcs = solver.pcs(npcs=2, pcscaling=1)
variance_fractions = solver.varianceFraction()
north_test = solver.northTest(vfscaled=True)


# --- spatial patterns
fig, ax = plot.subplots(axwidth=5, nrows=2, tight=True, proj='pcarree',
                        proj_kw={'lon_0': 180})
# --- format options
ax.format(land=False, coast=True, innerborders=True, borders=True,
          large='15px', labels=False,
          latlim=(31, -31), lonlim=(119, 291),
          geogridlinewidth=0,
          abcloc='ul')
# a) first EOF mode
map1 = ax[0].contourf(dset['lon'], dset['lat'], eofs[0, :, :],
                      levels=np.arange(-0.5, 0.6, 0.1), cmap='Div', extend='both')

# b) second EOF mode
map2 = ax[1].contourf(dset['lon'], dset['lat'], eofs[1, :, :],
                      levels=np.arange(-0.5, 0.6, 0.1), cmap='Div', extend='both')

ax[1].colorbar(map2, loc='b', length=0.9)

plt.show()


# --- principal components
fig, ax = plot.subplots(figsize=(10, 5), nrows=2, ncols=1,
                        tight=True, sharex=False, sharey=False)

ax[0].plot(np.arange(len(variance_fractions[0:10])),
           variance_fractions[0:10] * 100,)
ax[0].plot(np.arange(len(variance_fractions[0:10])),
           variance_fractions[0:10] * 100, 'ro')

ax[1].plot(pcs['time'], pcs[:, 0], color='red', label='PC1')
ax[1].plot(pcs['time'], pcs[:, 1], color='blue', label='PC2')

plt.legend()

plt.show()
