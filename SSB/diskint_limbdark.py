# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 14:01:19 2018

@author: foshe
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import functions as func

h, tau5, colm, temp, vturb, nhyd, nprot, nel, ptot, pgasptot,dens = np.loadtxt('falc.dat', usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)
wav, F, Fcont, I, Icont = np.loadtxt('solspect.dat', usecols=(0,1,2,3,4), unpack = True)

#Calculating intensity for the wavelengths given by table 5 in the assignment
int_calc = np.zeros(len(wav))
for i in range(len(wav)):
    hmean,int_calc[i],contfunc, tau = func.emergentIntensity(wav[i],h,tau5,temp, nel, nhyd, nprot,1)
    

#Comparing the calculated intensity given the FALC data with the observed values of Allen (1976)  
plt.figure(0)
plt.grid(linestyle='dotted')
plt.plot(wav,Icont, label ='Observed (Allen 1976)')
plt.plot(wav,int_calc*1e-14, label = 'Computed from FALC')
plt.xlabel('Wavelength [$\mu$m]', size = 12)
plt.ylabel('Intensity [$10^{14}$erg cm$^{-2}$ s$^{-1}$ cm$^{-1}$ ster$^{-1}$]', size = 12)
plt.title('Observed and computed solar disk-center intensity')
plt.xticks(np.arange(0, 5, 0.5))
plt.yticks(np.arange(0,6,0.5))
plt.legend()


#Computing intensity from different viewing angles
plt.figure(1)
plt.grid(linestyle='dotted')
plt.xticks(np.arange(0, 5, 0.5))
plt.yticks(np.arange(0,6,0.5))
plt.xlabel('Wavelength [$\mu$m]', size = 12)
plt.ylabel('Intensity [$10^{14}$erg cm$^{-2}$ s$^{-1}$ cm$^{-1}$ ster$^{-1}$]', size = 12)
plt.title('Computed emergent intensity from different viewing angles')
for mu in range(10):
    int_calc = np.zeros(len(wav))
    for i in range(len(wav)):
        hmean,int_calc[i],contfunc, tau = func.emergentIntensity(wav[i],h,tau5,temp, nel, nhyd, nprot, 0.1*(mu+1))
    plt.plot(wav,int_calc*1e-14, label = '$\mu = $' + '%.1f'% (0.1*(mu+1)))
plt.legend()


#Finding limb darknening for different wavelengths
wavelengths = [0.5,1,1.6,5]
mu = np.arange(0.01,1.01,0.01)
sintheta=np.sqrt(1-mu**2)

int_calc = np.zeros(len(mu))

plt.figure(2)
plt.grid(linestyle='dotted')
plt.xlabel('$\mu$', size = 12)
plt.ylabel('Normalized intensity', size = 12)
plt.title('Limb darkening for different wavelengths', size = 12)
plt.xticks(np.arange(0,1.1,0.1))
plt.yticks(np.arange(0,1.1,0.1))
plt.xlim(1,0)
plt.ylim(0,1)

for i in range(len(wavelengths)):
    for j in range(len(mu)):
        hmean,int_calc1,contfunc, tau = func.emergentIntensity(wavelengths[i],h,tau5,temp, nel, nhyd, nprot, 1)
        hmean,int_calc[j],contfunc, tau = func.emergentIntensity(wavelengths[i],h,tau5,temp, nel, nhyd, nprot, mu[j])
    plt.plot(mu,int_calc/int_calc1,label = '$\lambda= $' + '%.1f'% wavelengths[i] + ' $\mu$m')
#   plt.gca().invert_xaxis()


plt.legend()    

plt.figure(3)
plt.grid(linestyle='dotted')
plt.xlabel(r'$r/R_{\odot}$', size = 12)
plt.ylabel('Normalized intensity', size = 12)
plt.title('Limb darkening for different wavelengths', size = 12)
plt.xticks(np.arange(0,1.1,0.1))
plt.yticks(np.arange(0,1.1,0.1))
plt.xlim(0,1)
plt.ylim(0,1)
for i in range(len(wavelengths)):
    for j in range(len(mu)):
        hmean,int_calc1,contfunc, tau = func.emergentIntensity(wavelengths[i],h,tau5,temp, nel, nhyd, nprot, 1)
        hmean,int_calc[j],contfunc, tau = func.emergentIntensity(wavelengths[i],h,tau5,temp, nel, nhyd, nprot, mu[j])
    plt.plot(sintheta,int_calc/int_calc1,label = '$\lambda= $' + '%.1f'% wavelengths[i] + ' $\mu$m')


plt.legend()  
    
    