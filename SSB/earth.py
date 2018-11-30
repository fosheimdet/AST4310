# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 22:32:03 2018

@author: foshe
"""

import numpy as np
import matplotlib.pyplot as plt

#rc('font',**{'family':'serif'}) # This is for Latex writing
# DEFINE ALL THE CONSTANTS YOU NEED HERE (or any another place,
# if you prefer)

k=1.38065e-16    #Boltzmann constant
g_E = 980.665    #Gravitational constant on Earth
mH = 1.67352e-24 #Hydrogen mass
mHe = 3.97*mH    #Helium mass
mP = 1.67262e-24 #Proton mass

h, logP, temp, log_rho, logN = np.loadtxt('earth.dat',usecols=(0,1,2,3,4), unpack = True)


pressure = 10**logP
density = 10**log_rho
N = 10**logN

plt.figure(1)
plt.grid(linestyle='dotted')
plt.plot(h,temp)
plt.xlabel('Height [km]', size = 12)
plt.ylabel('Temperature [K]', size = 12)
plt.title('Temperature stratification of the Earth\'s atmosphere')

plt.figure(2)
plt.grid(linestyle='dotted')
plt.semilogy(h,10**(logP))
plt.xlabel('Height [km]', size = 12)
plt.ylabel('Pressure [dyn cm $^{-2}$]', size = 12)
plt.title('Pressure stratification of Earth\'s atmosphere')

plt.figure(3)
plt.grid(linestyle='dotted')
plt.semilogy(h,10**(log_rho))
plt.xlabel('Height [km]', size = 12)
plt.ylabel('Gas density [g cm $^{-2}$]', size = 12)
plt.title('Gas density in Earth\'s atmosphere')

plt.figure(4)
plt.grid(linestyle='dotted')
plt.semilogy(h,10**(logN))
plt.xlabel('Height [km]', size = 12)
plt.ylabel('Number density [cm $^{-3}$]', size = 12)
plt.title('Number density in Earth\'s atmosphere')


#Normalizing pressure and density
norm_P = pressure/np.max(pressure)
norm_rho=density/np.max(density)

plt.figure(5)
plt.grid(linestyle='dotted')
plt.plot(h,norm_P, label = 'Pressure')
plt.plot(h,norm_rho, label = 'Density')
plt.xlabel('Height [km]', size = 12)
plt.ylabel('Normalized pressure and density', size = 12)
plt.title('Pressure and density stratifications of Earth\'s atmosphere')

#Mean molecular weight
mu_E = density/(N*mH)
plt.figure(6)
plt.grid(linestyle='dotted')
plt.plot(h,mu_E)
plt.xlabel('Height [km]', size = 12)
plt.ylabel('Mean molecular weight in m$_H$', size = 12)
plt.title('Mean molecular weight of Earth\'s atmosphere')

#Finding scale height H_rho
#h[0]=0km
H_rho_ideal =0   #Scale height if the atmpshere were an ideal gas
#H_rho_exp = -h[1]/(np.log(density[1]/density[0]))  #Scale height determined from data
#H_rho = -h[1]/(np.log(density[1]/density[0]))
H_rho_exp=0
for i in range (1,9):
    print(i)
    H_rho_exp += -h[i]/(np.log(density[i]/density[0])*(8))
    H_rho_ideal+=k*temp[i]/((mu_E[i]*mH *g_E)*8)

#H_rho = H_rho/len(h)

print('Experimental scale height [km]: ', H_rho_exp, 'Ideal gas scale height [km]: ', H_rho_ideal/(1e2*1e3))

#Parameter values at the base of the atmosphere
print('Different values at the base of the atmosphere: ')
print('Pressure: ', pressure[0])
print('Temperature: ', temp[0])
print('Density: ', density[0])
print('Number density: ', N[0])

print('Column mass of Earth\'s atmosphere: ', pressure[0]/g_E)
print('Thermal photon production in the air around us: ', (20*temp[0]**3))

