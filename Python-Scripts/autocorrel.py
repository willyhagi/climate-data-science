#!/usr/bin/env python


'''
 File Name: autocorrel.py
 Description: Autocorrelation Function and Correlogram plot.
 Observations: The pandas.core.datetools module is deprecated, but
 this won't give you a problem for the calculations here.
 Author: Willy Hagi
 E-mail: hagi.willy@gmail.com
 Python Version: 3.6
'''





import numpy                    as np
import matplotlib.pylab         as plt

from statsmodels.tsa.stattools  import acf
from scipy.signal               import detrend



def dof(acfunc, confint, N):
    # estimating the effective sample size of data with acf
    lag  =  np.where(acfunc < conf[:,1])
    lag  =  lag[0][0]
    nef  =  N / lag
    return nef





##----- reading and detrending (optional) data 
sst  =  np.loadtxt('nino3.asc')
#sst  =  detrend(sst, type='linear')  # linear detrending



##----- autocorrelation function and confidence intervals
acf_total, conf  =  acf(sst, alpha=0.05, nlags=425, fft=True,
                        unbiased=True)
conf -= conf.mean(1)[:,None]  # lag 0 confidence intervals



##----- estimating the effective sample size
n   = len(sst) / 12. # number of years
print (n)
df  =  dof(acf_total, conf, n )
print (df)





##----------------------- PLOTTING

##--- CORRELOGRAM
plt.stem(np.arange(0,10,1), acf_total[0:10], linefmt='b-',
                            markerfmt='bo', basefmt='r-')
plt.plot(np.arange(0,10,1), conf[0:10], 'k--')
plt.title(u'Correlogram for Niño 3 Index')
plt.xlabel(u'Lags')
plt.ylabel(u'ACF')
plt.show()





# A. M .D .G #
