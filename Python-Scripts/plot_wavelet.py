#!/usr/bin/env python


import numpy                 as np
import matplotlib.pylab      as plt
import matplotlib.gridspec   as gridspec
import math
import matplotlib

from waveletFunctions        import*
from mpl_toolkits.axes_grid1 import make_axes_locatable


__author__ = 'Evgeniya Predybaylo'


'''
Minor modifications for Python 3.6 compatibility made by:
Willy Hagi (UEA/EST)
hagi.willy@gmail.com
'''

# WAVETEST Example Python script for WAVELET, using NINO3 SST dataset
#
# See "http://paos.colorado.edu/research/wavelets/"
# The Matlab code written January 1998 by C. Torrence is modified to Python by Evgeniya Predybaylo, December 2014
#
# Modified Oct 1999, changed Global Wavelet Spectrum (GWS) to be sideways,
#   changed all "log" to "log2", changed logarithmic axis on GWS to
#   a normal axis.
# ------------------------------------------------------------------------------------------------------------------



##--- autocorrelation function
def autocorrelation(x):
	'''
	from url: https://stackoverflow.com/questions/643699/how-can-i-use-numpy-correlate-to-do-autocorrelation
	'''
	xp = x-np.mean(x)
	f = np.fft.fft(xp)
	p = np.array([np.real(v)**2+np.imag(v)**2 for v in f])
	pi = np.fft.ifft(p)
	return np.real(pi)[:int(x.size/2)]/np.sum(xp**2)





##----- read data
sst_total = np.loadtxt('sst_nino3.dat')





#----------C-O-M-P-U-T-A-T-I-O-N------S-T-A-R-T-S------H-E-R-E------------------------------------------------------
'''normalize by standard deviation (not necessary, but makes it easier
   to compare with plot on Interactive Wavelet page, at
   http://paos.colorado.edu/research/wavelets/plot/'''

#variance  =  np.std(sst_total,ddof=0)**2
#sst = (sst - np.mean(sst)) / np.std(sst, ddof=1)
variance  =  np.std(sst_total, ddof=1) ** 2
stdev     =  np.std(sst_total,ddof=0)   # ddof = delta degrees of freedom. default is 0
sst       =  sst_total/ stdev           # normalized time series



##--- lag-1 autocorrelation for red noise background
acf     =  autocorrelation(sst)
alpha1  =  np.around(acf[1],decimals=2)
alpha2  =  np.around(acf[2],decimals=2)
lagf    =  (alpha1+math.sqrt(alpha2))/2.
lag1    =  np.around(lagf,decimals=2)    # lag-1 autocorrelation for red noise background



n    =  len(sst)
dt   =  (1/4)
time =  (np.arange(n))*dt + 1871.0   # construct time array
xlim =  ([1871, 2000])               # plotting range
pad  =  1.                           # pad the time series with zeroes (recommended)
dj   =  0.25                         # this will do 4 sub-octaves per octave
s0   =  2 * dt                       # this says start at a scale of 6 months
j1   =  7 / dj                       # this says do 7 powers-of-two with dj sub-octaves each
mother = 'MORLET'                    # options: 'MORLET', 'PAUL' e 'DOG'



##--- Wavelet Transform:
wave, period, scale, coi  = wavelet(sst, dt, pad, dj, s0, j1, mother)
power = (np.abs(wave))  ** 2.        # compute wavelet power spectrum



##---  Significance levels: (variance=1 for the normalized SST)
signif = wave_signif(([1.00]), dt=dt, sigtest=0, scale=scale, lag1=lag1, mother=mother)
sig95 = signif[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])  # expand signif --> (J+1)x(N) array
sig95 = power / sig95  # where ratio > 1, power is significant


##--- Global wavelet spectrum & significance levels:
global_ws = 1.00 * (np.sum(power, axis=1) / n)  # time-average over all times
dof = n - scale  # the -scale corrects for padding at edges
global_signif = wave_signif(1.00, dt=dt, scale=scale, sigtest=1, lag1=lag1, dof=dof, mother=mother)



##--- Scale-average between El Nino periods of 2--8 years
avg         =  np.logical_and(scale >=2., scale < 8.)
Cdelta      =  0.776  # this is for the MORLET wavelet
scale_avg   =  scale[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])  # expand scale --> (J+1)x(N) array
scale_avg   =  power / scale_avg  # [Eqn(24)]
scale_avg   =  1.00 * dj * dt / Cdelta * sum(scale_avg[avg, :])  # [Eqn(24)]
scaleavg_signif = wave_signif(1.00, dt=dt, scale=scale, sigtest=2, lag1=lag1, dof=([2, 8.]), mother=mother)



##--- RECTIFICATION OF BIAS (Liu et al, 2007)
powers      =   np.zeros_like(power)
for k in range(len(scale)):
    powers[k,:] = power[k,:] / scale[k]


global_wss      =   global_ws / scale
global_signifs  =   global_signif  / scale





#------------------------------------------------------ Plotting
gs = gridspec.GridSpec(2, 3)

plt.figure(figsize=(12,7))
##----- a) Contour plot wavelet power spectrum
plt3  = plt.subplot(gs[0,0:2])
levels = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16]
CS = plt.contourf(time, period, power, len(levels))  #*** or use 'contour'
#im = plt.contourf(CS, cmap='binary',vmin=1.0,vmax=8.) # more contour options
im  =  plt.contourf(CS, cmap='rainbow')
plt.xlabel('Time (years)',fontsize=11)
plt.ylabel('Period (years)',fontsize=11)
plt.title('a) Wavelet Power Spectrum',fontsize=11)
plt.xticks(range(1871, 2005, 25))

# 95% significance contour, levels at -99 (fake) and 1 (95% signif)
plt.contour(time, period, sig95, [-99, 1], colors='k')
# cone-of-influence, anything "below" is dubious
plt.plot(time, coi, 'k')

# 'x' fill under cone of influence area
ts = time
coi_area = np.concatenate([[np.max(scale)], coi, [np.max(scale)],[np.max(scale)]])
ts_area = np.concatenate([[ts[0]], ts, [ts[-1]] ,[ts[0]]])
plt.fill(ts_area,(coi_area),color='Black',alpha=0.5,hatch="x")

# format y-scale
plt3.set_yscale('log', basey=2, subsy=None)
plt.ylim([np.min(period), np.max(period)])
ax = plt.gca().yaxis
ax.set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt3.ticklabel_format(axis='y', style='plain')
plt3.invert_yaxis()

# set up the size and location of the colorbar
divider = make_axes_locatable(plt3)
cax = divider.append_axes("bottom", size="5%", pad=0.5)
plt.colorbar(im, cax=cax, orientation='horizontal')



##--- b) global wavelet spectrum
plt4  = plt.subplot(gs[0,2])
plt.plot(global_ws, period, 'black')
plt.plot(global_signif,period, 'k--')
plt.xlabel('Power (degC^2)',fontsize=11)
plt.ylabel('Periodo (years)',fontsize=11)
plt.title('b) Global Wavelet Spectrum',fontsize=11)
plt.xlim([0, 1.5 * np.max(global_ws)])
plt4.set_yscale('log', basey=2, subsy=None)
pmin = np.min(period)
pmax = np.max(period)
plt.ylim([np.min(period), np.max(period)])
ax = plt.gca().yaxis
ax.set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt4.ticklabel_format(axis='y', style='plain')
plt4.invert_yaxis()



##--- c) 2--8 yr scale-average time series
plt.subplot(gs[1,:])
plt.plot(time, scale_avg, 'black')
plt.xticks(range(1871, 2005, 25))
plt.xlabel('Time (years)',fontsize=11)
plt.ylabel('Avg variance (degC^2)',fontsize=11)
plt.title('c) 2-8 yr Scale-average Time Series',fontsize=11)
plt.plot(xlim, scaleavg_signif + [0, 0], 'k--')



plt.tight_layout()
plt.show()




# end of code
