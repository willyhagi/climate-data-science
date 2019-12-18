#!/usr/bin/env python


'''
 File Name: prime_meridian.py
 Description: Functions for operations with the prime meridian line.
 Author: Willy Hagi
 E-mail: hagi.willy@gmail.com
 Python Version: 3.7.3
'''


import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import cartopy as cartopy
import numpy as np
import xarray as xr
import proplot as plot

from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.util import add_cyclic_point


def fill_meridian(dataset, variable):
       '''
       Function to complete a gridded data without the prime meridian line.
       Currently, this is one of the most upsetting problems with xarray/cartopy/etc.

       Arguments:
       dataset: a gridded dataset.
       variable: your variable string.

       Returns:
       dset: a new gridded dataset with data in the prime meridian.
       '''
    lon_idx = dataset[variable].dims.index('lon')
    dset_c, lon = add_cyclic_point(dataset[variable].values,
                                   coord=dataset['lon'], axis=lon_idx)
    dset = xr.Dataset(data_vars={'lat': ('lat', dataset['lat']),
                                 'lon': ('lon', lon),
                                 'time': ('time', dataset['time']),
                                 variable: (['time', 'lat', 'lon'], dset_c)
                                 })
    return dset


def roll_meridian(dataset, central_lon, lon1, lon2):
    '''
    Function to get a longitude slice across the prime meridian line.
    Again, this is one of the major nuisances you can find.

    Arguments:
    dataset: a gridded dataset.
    central_lon: your central longitude (important when slicing).
    lon1: your selected first longitude of the new dataset.
    lon2: the same, but for the final longitude (across the prime meridian).

    Returns:
    dset: a new sliced dataset.
    '''
    lon = dataset['lon']
    lonroll = lon.roll(lon=central_lon)
    lonroll = lonroll.sel(lon=slice(lon1, lon2))
    dset = dataset.sel(lon=lonroll)
    return dset

# --- read netcdf file
dset = xr.open_dataset('sst.mnmean.nc').load()

# --- use function to complete the dataset
prime_fill = fill_meridian(dset, 'sst')

# --- use function to slice across the prime meridian longitude
prime_slice = roll_meridian(prime_fill['sst'], 150, 100., 20.)


# --- plot the comparison between the two datasets
fig, ax = plot.subplots(axwidth=5, nrows=2, ncols=1,
                        tight=True, proj='pcarree')

ax.format(land=False, coast=True, innerborders=True, borders=True, grid=False, labels=True)

ax[0].contourf(dset['lon'], dset['lat'], dset['sst'][0, :, :], cmap='Marine')
ax[1].contourf(prime_fill['lon'], prime_fill['lat'],
               prime_fill['sst'][0, :, :], cmap='Marine')


plt.show()


# --- plot the comparison between the two datasets (again)
fig, ax = plot.subplots(axwidth=5, nrows=2, ncols=1,
                        tight=True, proj='pcarree', proj_kw={'lon_0': 0})
ax.format(land=False, coast=True, innerborders=True, borders=True, grid=False, labels=True)

ax[0].contourf(prime_fill['lon'], prime_fill['lat'],
               prime_fill['sst'][0, :, :], cmap='Marine')
ax[1].contourf(prime_slice['lon'], prime_slice['lat'], prime_slice[0,:,:], cmap='Marine')
plt.show()
