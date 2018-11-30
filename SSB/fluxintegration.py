# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 19:50:28 2018

@author: foshe
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import functions as func


h, tau5, colm, temp, vturb, nhyd, nprot, nel, ptot, pgasptot,dens = np.loadtxt('falc.dat', usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)
wav, F, Fcont, I, Icont = np.loadtxt('solspect.dat', usecols=(0,1,2,3,4), unpack = True)

sigma_Thomson=6.648e-25
# SSB 2.7 page 17: flux integration
# ===== three-point Gaussian integration intensity -> flux
# abscissae + weights n=3 Abramowitz & Stegun page 916
xgauss=[-0.7745966692,0.0000000000,0.7745966692]
wgauss=[ 0.5555555555,0.8888888888,0.5555555555]
fluxspec = np.zeros(len(wav),dtype=float)
intmu = np.zeros((3,len(wav)), dtype=float)
for imu in range(3):
    mu=0.5+xgauss[imu]/2. # rescale xrange [-1,+1] to [0,1]
    wg=wgauss[imu]/2. # weights add up to 2 on [-1,+1]
    for iw in range(0,len(wav)):
        wl=wav[iw]
        ext = np.zeros(len(tau5))
        tau = np.zeros(len(tau5))
        integrand = np.zeros(len(tau5))
        intt = 0.0
        for i in range(1, len(tau5)):
            ext[i] = (func.exthmin(wl*1e4, temp[i], nel[i])*(nhyd[i]-nprot[i]) + sigma_Thomson*nel[i])
            tau[i] = (tau[i-1] + 0.5 * (ext[i] + ext[i-1]) *
            (h[i-1]-h[i])*1E5)
            integrand[i] = func.planck(temp[i],wl*1e-4)*np.exp(-tau[i]/mu)
            intt += 0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])/mu
        intmu[imu,iw]=intt
        fluxspec[iw]=fluxspec[iw] + wg*intmu[imu,iw]*mu
fluxspec *= 2 # no np.pi, Allen 1978 has flux F, not {\cal F}
figname='ssb_2.7_fluxintegration'
f=plt.figure(figname)
plt.plot(wav,fluxspec*1e-14, label='computed from FALC')
plt.plot(wav,Fcont, label='observed (Allen 1978)')
plt.legend(loc='upper right')
plt.title('observed and computed continuum flux')
plt.ylabel(r'astrophysical flux [10$^{14}$ erg s$^{-1}$ cm$^{-2}$ ster$^{-1}$ cm$^{-1}$]')
plt.xlabel('wavelength [$\mu$m]')
plt.xticks(np.arange(0,5,0.5))
plt.grid(linestyle='dotted')
plt.show()
f.savefig(figname+'.pdf',bbox_inches='tight')
f.savefig(figname+'.png',bbox_inches='tight')