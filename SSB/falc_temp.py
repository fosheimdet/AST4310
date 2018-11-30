# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 13:37:07 2018

@author: foshe
"""

import numpy as np # numerical package
import matplotlib.pyplot as plt # plotting package
from matplotlib import rc
rc('font',**{'family':'serif'}) # This is for Latex writing
# DEFINE ALL THE CONSTANTS YOU NEED HERE (or any another place,
# if you prefer)
# reading falc.dat
k=1.38065e-16    #Boltzmann constant
mH = 1.67352e-24 #Hydrogen mass
mHe = 3.97*mH    #Helium mass
mP = 1.67262e-24 #Proton mass



h, tau5, colm, temp, vturb, nhyd, nprot, nel, ptot, pgasptot,dens = np.loadtxt('falc.dat', usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)

#(h, tau5, colm, temp, vturb, nhyd, nprot, nel, ptot, pgasptot,
#dens = np.loadtxt('falc.dat', usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)) 
# plotting



fig0 = plt.figure(0)
plt.title('Temperature stratification of the Sun\'s atmosphere')
plt.grid(linestyle='dotted')
plt.plot(h, temp)
plt.xlabel('height [km]', size = 12)
plt.ylabel('temperature [K]', size = 12)
#fig.savefig('Myfigure.pdf', bbox_inches='tight',
#pad_inches=0.106)
plt.xlim(-500,2500)
plt.ylim(2000,10000)


#Finding the gravitational constant
g_S = ptot[55]/colm[55]
print('The gravitational constant on the surface is: ', g_S)

fig1 = plt.figure(1)
plt.title('Pressure plotted against column mass ')
ax1 = fig1.add_subplot(2,1,1)
ax2 = fig1.add_subplot(2,1,2)
ax1.plot(colm, ptot)
ax2.plot(colm,np.log(ptot))
ax1.set_ylabel('total pressure [dyn $\mathrm{cm}^{-2}$]', size = 12)
#ax1.set_xlabel('column mass [g $\mathrm{cm}^{-2}$]',size = 12)

ax2.set_ylabel('log(total pressure) [dyn $\mathrm{cm}^{-2}$]',size = 12)
ax2.set_xlabel('column mass [g $\mathrm{cm}^{-2}$]',size = 12)

#Finding the mass densities of Hydrogen and Helium as well the density fraction of the metals
rho_H = nhyd*mH
rho_He= nhyd*0.1*mHe
rho_remaining = 1 - np.average((rho_H + rho_He)/dens)
print('Mass density fraction of metals: ',rho_remaining)


plt.figure(2)
plt.title('Hydrogen and Helium abundance')
plt.grid(linestyle='dotted')
plt.plot(h, rho_H/dens, label ='Hydrogen')
plt.plot(h, rho_He/dens, label = 'Helium')
plt.xlabel('height [km]', size = 12)
plt.ylabel('mass density fraction', size = 12)
plt.legend()

#Plotting column mass against height
fig3 = plt.figure(3)
plt.grid(linestyle='dotted')
ax1 = fig3.add_subplot(2,1,1)
ax2 = fig3.add_subplot(2,1,2)
ax1.plot(h,colm)
ax2.plot(h,np.log(colm))
ax1.set_ylabel('column mass [g cm$^{-2}$]')

ax2.set_ylabel('log(column mass) [g cm$^{-2}$]')
ax2.set_xlabel('height [km]')

#Find index of element in h with h=0

hindex = 0
for i in range(len(h)):
    if(h[i]==0):
        hindex = i

#Estimating the scale height in the deep photosphere
H_rho_ideal = (k*temp[-1]/(mH*g_S))/(1e2*1e3)

H_rho_exp =-(h[20]/np.log(dens[20]/dens[hindex])) 
#
#print(h) 

print('Density scale height, ideal gas: ', H_rho_ideal, 'Density scale height, experimental: ', H_rho_exp)
#H_rho = 146.2
#print(H_rho)
plt.figure(4)
plt.grid(linestyle='dotted')
plt.plot(h,dens, label = 'FALC')
plt.plot(h,dens[hindex]*np.exp(-h/H_rho_ideal), label = r'H$_{{\rho}_{ideal}} = 283.1 km$')
plt.plot(h,dens[hindex]*np.exp(-h/H_rho_exp), label = r'H$_{{\rho}_{exp}} = 147.1 km$')
plt.xlabel('Height [km]', size = 12)
plt.ylabel('Gas density [g cm $^{-3}$]', size = 12)
plt.title('FALC vs Ideal gas')
plt.legend()

Pfalc = pgasptot*ptot   #FALC gas pressure
Pideal = (nhyd+nel)*k*temp
Pideal_he = (nhyd+nel+0.1*nhyd)*k*temp

plt.figure(5)
plt.grid(linestyle='dotted')
plt.plot(h,(Pfalc), label = 'FALC')
plt.plot(h,(Pideal), label = '$(n_H + n_e)kT$')
plt.title('Gas pressure against height')
plt.xlabel('Height [km]')
plt.ylabel('Gas pressure [dyne cm $^{-2}$]')
plt.legend()

fig6 = plt.figure(6)
ax1 = fig6.add_subplot(2,1,1)
ax1.set_title('Gas pressure ratios, Helium excluded')
ax2 = fig6.add_subplot(2,1,2)
ax2.set_title('Gas pressure ratios, Helium included')
ax1.plot(h,Pfalc/Pideal)
ax2.plot(h,Pfalc/Pideal_he)
ax1.set_ylabel('Gas pressure ratio', size = 12)
ax2.set_ylabel('Gas pressure ratio', size = 12)
ax2.set_xlabel('height [km]', size = 12)

#Various number densities plotted vs height
plt.figure(7)
plt.grid(linestyle='dotted')
plt.semilogy(h,(nhyd),label = "Hydrogen density")
plt.semilogy(h,(nel), label = "Electron density")
plt.semilogy(h,nprot, label = "Proton density")
plt.semilogy(h,(nel-nprot),label ="Electrons not from H")
plt.ylabel('Number density [cm$^{-3}$]', size = 12)
plt.xlabel('Height [km]', size = 12)
plt.title('Various number densities plotted vs height')
plt.legend()

#Hydrogen ionization fraction
H_ionfrac = (nprot)/(nprot+nhyd)
plt.figure(8)
plt.grid(linestyle='dotted')
plt.title('Hydrogen ionization fraction')
plt.semilogy(h,H_ionfrac)
plt.xlabel('Height [km]', size = 12)
plt.ylabel('Ionization fraction of H', size = 12)

#Photon density vs hydrogen density
#Valid for the deepest model location (h=-100km)

N_photon_lowest= 20*temp[-1]**3
print('Photon density at -100km: ', N_photon_lowest, 'Hydrogen density at -100km: ', nhyd[-1])

N_photon_highest = 20*(5770**3)/(2*np.pi)
print('Photon density at 2218km: ', N_photon_highest, 'Hydrogen density at 2218km: ', nhyd[0])

#Particle density at h=0
print('Particle density at h=0: ', (nhyd[hindex]+nel[hindex]+0.1*nhyd[hindex]))

print('Column mass at h=0, sun: ', ptot[hindex]/g_S)

