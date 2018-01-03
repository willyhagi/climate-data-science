#!/usr/bin/env python


'''
 File Name: rad_laws.py
 Description: Introduction to Classes with Radiation Laws.
 Observations: None.
 Author: Willy Hagi
 E-mail: hagi.willy@gmail.com
 Python Version: 3.6
'''





import numpy as np
import matplotlib.pyplot as plt
from scipy.special import expit





##----- DEFINING THE CLASS
class Rad:
    def __init__(self, T, wvlg, epsilon): 
        self.T        =  T            # temperature (K)
        self.wvlg     =  wvlg * 1E-6  # wavelength (m) 
        self.epsilon  =  epsilon      # emissivity
       
        
    def planck(self, c1 = 3.74E-16, c2 = 1.45E-2):
        '''
        Planck Function

        [c1]  =   W m^2      
        [c2]  =   m K
        First and second radiative constants, respectively.
        '''
        
        self.c    =  c2 / (self.wvlg * self.T) # this is just for simplicity
        self.ex   =  np.exp(self.c, dtype=np.float64)
        self.B    =  c1 * (self.wvlg ** -5.) / (np.pi*(self.ex-1))
        return  self.B # intensity of emitted radiation

    def stefan(self, sigma = 5.67E-8):
        '''
        Stefan-Boltzmann Law

        sigma   =  stefan-boltzmann constant (W m^-2 K^-4)
        '''
        self.E  =  self.epsilon * sigma * (self.T**4.)
        return self.E

    def wien(self):
        '''
        Wien displacement Law
        '''

        self.wmax = 2897 / self.T  # wavelength (micrometer) of the maximum monocromatic intensity
        return self.wmax
        
    


     
lmbd = np.arange(0.1,2.0,0.01)  # wavelength (micrometer)


## class arguments
B1  =  Rad(5000,lmbd,1) 
B2  =  Rad(6000,lmbd,1)  
B3  =  Rad(7000,lmbd,1) 


## planck's law
P1  =  B1.planck()
P2  =  B2.planck()
P3  =  B3.planck()


## stefan-boltzmann law
S1  =  B1.stefan()
S2  =  B2.stefan()
S3  =  B3.stefan()


## wien's law
W1  =  B1.wien()
W2  =  B2.wien()
W3  =  B3.wien()




##----- PLOTTING
plt.plot(lmbd, P1 * 1E-6,'blue',   label='5000 K')
plt.plot(lmbd, P2 * 1E-6,'orange', label='6000 K')
plt.plot(lmbd, P3 * 1E-6,'red',    label='7000 K')
plt.legend(loc='upper right')
plt.ticklabel_format(style='sci', axis='y')
plt.title("Planck's Function")
plt.xlabel(r'Wavelength ($\mu m$)')
plt.ylabel(r'$ B_{\lambda} (W m^{2} sr^{-1} \mu m) \times 10^{7}  $')


plt.show()





# A. M. D. G #
