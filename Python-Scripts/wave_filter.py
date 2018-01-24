#!/usr/bin/env python



'''
 File Name: wave_filter.py
 Description: Wavelet Filtering method.
 Observations: Your input data must be a time series / the wavelet
 functions were provided by Evgenyia Predybaylo and are available at:
 http://paos.colorado.edu/research/wavelets/software.html  .
 Author: Willy Hagi
 E-mail: hagi.willy@gmail.com
 Python Version: 3.6
'''





import numpy             as np
import matplotlib.pyplot as plt
import pandas            as pd
import math
import matplotlib

from waveletFunctions    import*





##----- read time series
sst  =  np.loadtxt('sst_nino3.dat')



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
realpart  =  np.real(wave)




##--- filtering
cdelta    =  0.776                              # cdelta=0.776 for morlet wavelet; see table 2 (TC, 1998)
psi       =  math.pi**(-1/4)                    # for morlet wavelet; see table 2 (TC, 1998)
xnc       =  (dj*math.sqrt(dt))/(cdelta*psi)    # 'constant' part for Eq. 29 (TC, 1998)
avg       =  np.logical_and(scale>=2.,scale<8.) # filtering scales
scale_avg =  scale[:, np.newaxis].dot(np.ones(n)[np.newaxis, :]) # expand scale array
xnv       =  realpart/np.sqrt(scale_avg)        # 'variable' part for Eq. 29
xn        =  xnc*sum(xnv[avg,:])                # Eq. 29


##--- save file
#np.savetxt('nino3_2to8yrs.asc',xn,fmt='%.10f')





##----------------------- PLOTTING
plt.figure(figsize=(12,4))
plt.xlim([1871,2000])
plt.xticks(range(1871, 2005, 25))
plt.plot(time, sst,'black',lw=1,label=u"Niño 3")
plt.plot(time,xn,'red',lw=1,label="2-8 yrs")
plt.xlabel('Years')
plt.ylabel(u'ASST ($\degree$C)')
plt.title(u'2-8 yrs Niño 3 Wavelet Reconstruction')
plt.legend(loc='upper right')

plt.tight_layout
plt.show()





# A. M. D. G #
