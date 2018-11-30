# -*- coding: utf-8 -*-
"""
Created on Sun Nov 25 18:10:35 2018

@author: foshe
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import functions.py as func
#from solspect import exthmin

h, tau5, colm, temp, vturb, nhyd, nprot, nel, ptot, pgasptot, dens = np.loadtxt("falc.dat",usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)
wav, F,Fcont,I,Icont= np.loadtxt("solspect.dat",usecols=(0,1,2,3,4), unpack=True)



tau = np.zeros(len(tau5), dtype=float) # initializing tau array 
nneutH=nhyd-nprot
sigma_T= 6.648e-25
EXmin = func.exthmin(5000,temp,nel) 
Thomson = nel*sigma_T
ext = EXmin+Thomson

for i in range(1,len(tau)):
	tau[i] = tau[i-1] + 0.5*(ext[i]+ext[i-1])*(h[i-1]-h[i])*1e5
# index zero is not accounted for, so tau[0] = 0 because we have already initialized

fig = plt.figure()
plt.grid(linestyle='dotted')
plt.plot(h,tau5, label = 'tau5')
plt.plot(h,tau, label = 'tau')
plt.yscale('log')
plt.xlabel(r'Height [km]')
plt.ylabel(r'Optical depth[$\tau_{\lambda}$]')
plt.legend([r'$\tau_{500}$ in FALC',r'$\tau_{500}$ from H- and Thomson extinction'])
#
#tau = np.zeros(len(tau5), dtype=float) # initializing tau array
#ext = exthmin(500e-9, temp, nel)
#for i in range(1,len(tau)):
#    tau[i] = tau[i-1] + 0.5*(ext[i]+ext[i-1])*(h[i-1]-h[i])*1e5
## index zero is not accounted for, so tau[0] = 0 because we have already initialized
#plt.plot(h,tau5,'--', label = 'tau5')
#plt.plot(h,tau, label = 'tau')
#plt.yscale('log')
#plt.show()