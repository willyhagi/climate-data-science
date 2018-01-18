#!/usr/bin/env python


'''
 File Name: functions.py
 Description: A list of useful python functions and classes.
 Observations: You'll need this file to run the other scripts
from the climate-statistics repository.
 Author: Willy Hagi
 E-mail: hagi.willy@gmail.com
 Python Version: 3.6
'''





import numpy               as np
import xarray              as xr

from scipy.stats          import pearsonr





##----------------------- CLASSES


##----- monthly gridded averages, anomalies and indices

class Demean:
    def __init__(self, data, lat1, lat2, lon1, lon2, time1, time2):
        self.data  =  data
        self.lat1  =  lat1  ; self.lat2  = lat2
        self.lon1  =  lon1  ; self.lon2  = lon2
        self.time1 =  time1 ; self.time2 = time2


    def average(self):
        '''
        Monthly gridded mean values.
        Functions 'sel' and 'groupby' are from xarray package.
        '''
        self.x     =  self.data.sel(time  =  slice(self.time1, self.time2),
                                    lat   =  slice(self.lat1,  self.lat2),
                                    lon   =  slice(self.lon1,  self.lon2))
                               #lvl = lvl1
        self.xmean =  self.x.groupby('time.month').mean('time')
        return self.xmean

    def anomaly(self):
        # monthly gridded anomaly values

        self.anom = self.x.groupby('time.month') - self.xmean
        return self.anom

    def area_avg(self):
        '''
        This will make an average over space, if you wish to have an index.
        Please, notice that it will use as lat/lon values the ones you defined above
        for gridded calculations.
        '''

        self.index = self.anom.mean(dim = ('lat','lon'), skipna = True)
        self.index = np.asarray(self.index)
        return self.index



##----- select months for seasonal calculations

class Seasonal:
    def __init__(self, data, m1, m2, m3): # this can be extended to n months
        self.data  =  data
        self.m1    =  m1
        self.m2    =  m2
        self.m3    =  m3

    def select_months(self, month):
        return (month == self.m1) | (month == self.m2) | (month == self.m3)

    def season(self):
        # season will return n months for the whole time period of your dataset
        self.svar =  self.data.sel(time = self.select_months(self.data['time.month']))
        return self.svar



##----- masking NaN values from gridded data
class Nan_calc:
    def __init__(self, x, time, lat, lon):
        self.x     =  x
        self.time  =  time
        self.lat   =  lat
        self.lon   =  lon

    def reshaping(self):
        self.xnp   =  np.asarray(self.x)
        self.xres  =  np.reshape(self.xnp, (len(self.time), len(self.lat) * len(self.lon)), order='F')
        return self.xres

    def masking_array(self):
        self.masked_array  =  np.ma.masked_array(self.xres, np.isnan(self.xres))
        self.nan_values    =  self.masked_array.sum(0).mask
        self.values        =  ~self.nan_values
        self.x_masked      =  self.xres[:,self.values]
        return self.values, self.x_masked





##----------------------- FUNCTIONS

# gridded anomaly values
def anom(x,y):
    return x - y.sel(time=slice(time1,time2),
                     lat=slice(lat1,lat2),
                     lon=slice(lon1,lon2)).mean('time')


# matrix reconstruction
def rec_matrix(x, time, lat, lon, val):
    rec  =  np.ones( (len(time), len(lat)*len(lon)) ) * -999.
    for i in range(len(time)):
        rec[i,val]  =  x[i,:]
    rec_res    =  np.reshape(rec, (len(time), len(lat), len(lon)), order='F')
    rec_final  =  np.ma.masked_values(rec_res, -999.)
    return rec_final


# EOF using SVD
def eof_svd(F):
    U, g, V  =  np.linalg.svd(F, full_matrices=False)
    U        =  np.real(U)            # left singular vectors
    g        =  np.real(g)            # singular values
    V        =  np.real(V)            # right singular vectors
    #
    lmbd     =  g **2.                # eigenvalues
    PC       =  np.dot(np.diag(g), U) # principal components
    EOF      =  np.dot(V.transpose(), PC)         # EOF projection
    return U, V, g, lmbd, PC.transpose(), EOF


# EOF normalization
def eof_norm(PC, V, lmbd):
    PCnor   =  PC / np.sqrt(lmbd)
    EOFnor  =  V * np.sqrt(lmbd)
    return PCnor, EOFnor


# column-wise pearson linear correlation coefficient
def pearson(matrix, index, grid):
    cont = 0
    r    =  np.zeros(grid)
    p    =  np.zeros(grid)
    for i in np.arange(1, grid, 1):
        r[cont], p[cont] = pearsonr(matrix[:,i-1], index)
        cont += 1
    return r, p


#---- pearson correlation's student t-test
def t_test(r, dof, ttab):
    t_num   =  np.multiply(r, np.sqrt(dof - 2.) )
    t_den   =  np.sqrt(1. - r ** 2.)
    t_calc  =  np.divide(t_num, t_den)
    return t_calc / ttab

   
#--- composite analysis t-test
def t_comp(x, N, ttab):
    varx   =  np.sqrt(x.var() / N)
    tcalc  =  x / varx
    return tcalc / ttab  
