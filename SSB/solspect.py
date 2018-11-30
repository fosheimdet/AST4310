# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 21:55:45 2018

@author: foshe
"""

import numpy as np
import matplotlib.pyplot as plt
import functions as func



wav, F, Fcont, I, Icont = np.loadtxt('solspect.dat', usecols=(0,1,2,3,4), unpack = True)

c= 2.99792e10    #Speed of light in cm per second
h= 6.62607e-27   #Planck constant
k=1.38065e-16    #Boltzmann constant
g_E = 980.665    #Gravitational constant on Earth
mH = 1.67352e-24 #Hydrogen mass
mHe = 3.97*mH    #Helium mass
mP = 1.67262e-24 #Proton mass

factor = (wav**2/(c))*1e6#factor to convert intensity from wavelength domain to frequency domain

#def planck(temp, wav):
#    h = 6.62607e-27
#    k = 1.38065e-16
#    c = 2.99792e10
#    B = ((2*h*c**2)/(wav**5))*(1/(np.exp(h*c/(wav*k*temp))-1))
#    return B
#
#def BrightnessTemperature(wav, I):
#
#
#	#Inverts the Planck function and returns the brightness temperature
#
#	return ((h*c)/(k*wav))/np.log(1 + ((2*h*c**2)/(I*wav**5))*1e-4)


plt.figure(0)
#plt.grid("on")
plt.grid(linestyle='dotted')
plt.plot(wav,F, label = 'F$_{\lambda}$')
plt.plot(wav,Fcont, label = 'F\'$_{\lambda}$')
plt.plot(wav,I, label = 'I$_{\lambda}$')
plt.plot(wav,Icont, label = 'I\'$_{\lambda}$')
plt.xlabel('Wavelength [$\mu m$]', size = 12)
plt.ylabel('Intensity [$10^{10}$erg cm$^{-2}$ s$^{-1}$ $\mu m^{-1}$ ster$^{-1}$]', size = 12)
plt.xlim(0.2,5)
plt.xticks(np.arange(0, 5, 0.5))
plt.yticks(np.arange(0,5,0.5))
plt.title('Solar spectral distribution')
plt.legend()

print ('max(Ic)= ',np.max(Icont),'at',wav[np.where(Icont == np.max(Icont))])

#Solar spectrum per frequency bandwidth
plt.figure(1)
plt.grid(linestyle='dotted')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.plot(wav,F*factor, label = 'F$_{\lambda}$')
plt.plot(wav,Fcont*factor, label = 'F\'$_{\lambda}$')
plt.plot(wav,I*factor, label = 'I$_{\lambda}$')
plt.plot(wav,Icont*factor, label = 'I\'$_{\lambda}$')
plt.xlabel('Wavelength [$\mu m$]', size = 12)
plt.ylabel('Intensity [erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$ ster$^{-1}$]', size = 12)
plt.title('Solar spectral distribution ')
plt.legend()
plt.xticks(np.arange(0, 5, 0.5))
plt.yticks(np.arange(0,5e-5,0.5e-5))

print ('max(Ic)= ',np.max(Icont*factor),'at',wav[np.where(Icont*factor == np.max(Icont*factor))])

#Fitting the planck function to Icont
T_fit = 6330
plt.figure(2)
plt.grid(linestyle='dotted')
plt.plot(wav,func.planck(T_fit,wav/1e4)*1e-14, label ='Planck, T=6330')
plt.plot(wav,Icont, label = 'I\'$_{\lambda}$')
plt.xlabel('Wavelength [$\mu m$]', size = 12)
plt.ylabel('Intensity [$10^{10}$erg cm$^{-2}$ s$^{-1}$ $\mu m^{-1}$ ster$^{-1}$]', size = 12)
plt.xlim(0.2,5)
plt.xticks(np.arange(0, 5, 0.5))
plt.yticks(np.arange(0,5,0.5))
plt.title('Planck fit to the continuum intensity')
plt.legend()

#Finding brightness temperature
T_b = h*c/((wav*1e-4)*k*np.log(((2*h*c**2)/(((wav*1e-4)**5)*Icont*1e-10))+1))

plt.figure(3)
plt.grid(linestyle='dotted')
plt.plot(wav,func.BrightnessTemperature(wav*1e-4,Icont*1e10))
plt.xlabel('Wavelength [$\mu m$]', size = 12)
plt.ylabel('Brigthness temperature [K]')
plt.title('Brightness temperature as a function of wavelength')
plt.xticks(np.arange(0, 5, 0.5))
plt.yticks(np.arange(5200,6900,100))

