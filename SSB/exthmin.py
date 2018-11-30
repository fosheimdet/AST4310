# -*- coding: utf-8 -*-
"""
Created on Sun Nov 25 13:18:38 2018

@author: foshe
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import functions as func



h, tau5, colm, temp, vturb, nhyd, nprot, nel, ptot, pgasptot,dens = np.loadtxt('falc.dat', usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)
wav, F, Fcont, I, Icont = np.loadtxt('solspect.dat', usecols=(0,1,2,3,4), unpack = True)




#Plotting extinction at h=0
T = temp[np.argwhere(h ==0)[0][0]]
n_el = nel[np.argwhere(h ==0)[0][0]]
wavelength = np.linspace(wav[0],2,80)*1e4 #multiplying with 1e4 to convert from microns to ångstrøm
EXmin = func.exthmin(wavelength,T,n_el) 

plt.figure(0)
plt.title("H$^{-}$ extinction at h = 0")
plt.plot(wavelength*1e-4,EXmin)
plt.grid(linestyle='dotted')
plt.xlabel(r"Wavelength $\lambda$[$\mu$m]", size = 12)
plt.ylabel(r"H$^{-}$ extinction [cm$^2$ per H-atom]", size =12)

plt.figure(1)
plt.title("H$^{-}$ transmittance at h = 0")
plt.plot(wavelength*1e-4,1/EXmin)
plt.grid(linestyle='dotted')
plt.xlabel(r"Wavelength $\lambda$[$\mu$m]", size = 12)
plt.ylabel(r"H$^{-}$ transmittance [cm$^{-2}$ per H-atom]", size =12)

#plotting the extinction due to H- measure per cm path length
nneutH = nhyd-nprot
alpha_H = nneutH *func.exthmin(0.5e4,temp,nel)
plt.figure(2)
plt.grid(linestyle='dotted')
plt.semilogy(h,alpha_H)
plt.xticks(np.arange(-200, 2500, 250))
plt.xlabel('Height [km]', size = 12)
plt.ylabel(r'$\alpha_{\lambda}$ [cm$^{-1}]$]', size = 12)
plt.title(r'H$^-$ Extinction at $\lambda =0.5 \mu m$ as a function of height ')

#Adding Thomson scattering off free electrons

sigma_T = 6.648e-25
alpha_T = sigma_T*nel

plt.figure(3)
plt.grid(linestyle='dotted')
plt.semilogy(h,alpha_H, label = 'H$^{-}$')
plt.semilogy(h,alpha_T, label ='Thomson')
plt.semilogy(h,alpha_H+alpha_T, label = 'Total')
plt.xlabel('Height [km]', size = 12)
plt.ylabel(r'$\alpha_{\lambda}$ [cm$^{-1}]$]', size = 12)
plt.title(r'Extinction at $\lambda =0.5 \mu m$ as a function of height ')
plt.legend()



#Finding and plotting optical depth 

tau = np.zeros(len(tau5), dtype=float) # initializing tau array
alpha_H = func.exthmin(5000, temp, nel)*nneutH
sigma_T = 6.648e-25
alpha_T = sigma_T*nel
alpha_tot = alpha_H+alpha_T
for i in range(1,len(tau)):
    tau[i] = tau[i-1] + 0.5*(alpha_tot[i]+alpha_tot[i-1])*(h[i-1]-h[i])*1e5
# index zero is not accounted for, so tau[0] = 0 because we have already initialized
plt.figure(5)
plt.grid(linestyle='dotted')
plt.plot(h,tau5,'--', label = 'tau5')
plt.plot(h,tau, label = 'tau')
plt.yscale('log')
plt.xlabel('Height [km]')
plt.ylabel(r'Optical depth[$\tau_{\lambda}$]')
plt.legend([r'$\tau_{500}$ in FALC',r'$\tau_{500}$ from H- and Thomson extinction'])
plt.show()




hmean0,intt0,contfunc0, tau0 = func.emergentIntensity(0.5,h,tau5,temp, nel, nhyd, nprot, 1)
hmean1,intt1,contfunc1, tau1 = func.emergentIntensity(1,h,tau5,temp, nel, nhyd, nprot, 1)
hmean2,intt2,contfunc2, tau2 = func.emergentIntensity(1.6,h,tau5,temp, nel, nhyd, nprot, 1)
hmean3,intt3,contfunc3, tau3 = func.emergentIntensity(5,h,tau5,temp, nel, nhyd, nprot, 1)
print ('computed continuum intensity wl =%g : %g erg s-1 cm-2 ster-1 cm-1'
%(0.5, intt0))
w = np.where(wav == 0.5)
print ('observed continuum intensity wav=%g : %g erg s-1 cm-2 ster-1 cm-1'
%(wav[w], Icont[w]*1e10*1e4))

print('Mean height of formation', hmean0)


#peaknorm_cont = contfunc/np.max(contfunc)
plt.figure(6)
plt.grid(linestyle='dotted')
plt.plot(h,contfunc0/np.max(contfunc0), label = r'$\lambda$ = 0.5 $\mu$m')
plt.axvline(x=hmean0, color = "black",linestyle = "--", alpha = 0.7, label = r"$<h> = $%.1f km"%hmean0)
plt.xlabel('Height [km]')
plt.ylabel('Contribution function')
plt.title('Contribution function at $\lambda =0.5 \mu$m for different heights')
plt.xticks(np.arange(-100, 500, 100))
plt.xlim(-100,500)
plt.legend()
#plt.plot(hmean0*np.zeros(80),contfunc0/np.max(contfunc0))
#plt.plot([hmean0.value], [contfunc0[np.argmin(abs((h-hmean0).value))].value], marker='o',ms=9,color='blue')
#plt.plot(-3.1397,contfunc0[np.where(wav==-3.1397)]/np.max(contfunc0[np.where(wav==-3.1397)]) )

plt.figure(7)
plt.grid(linestyle='dotted')
plt.plot(h,contfunc0/np.max(contfunc0), label = r'$\lambda$ = 0.5 $\mu$m')
plt.plot(h,contfunc1/np.max(contfunc1), label = r'$\lambda$ = 1.0 $\mu$m')
plt.plot(h,contfunc2/np.max(contfunc2), label = r'$\lambda$ = 1.6 $\mu$m')
plt.plot(h,contfunc3/np.max(contfunc3), label = r'$\lambda$ = 5.0 $\mu$m')
plt.xticks(np.arange(-100, 500, 100))
plt.xlim(-100,500)
plt.xlabel('Height [km]', size = 12)
plt.ylabel('Contribution function', size = 12)
plt.title('Peak normalized contribution function')
plt.legend()

#tau = np.zeros(len(tau5), dtype=float) # initializing tau array 
#nneutH=nhyd-nprot
#sigma_T= 6.648e-25
#EXmin = exthmin(5000,temp,nel) 
#Thomson = nel*sigma_T
#ext = EXmin+Thomson
#
#for i in range(1,len(tau)):
#	tau[i] = tau[i-1] + 0.5*(ext[i]+ext[i-1])*(h[i-1]-h[i])*1e5
#print(tau)
print('500 nm: Mean height of formation: ', hmean0, 'h(tau=1): ', h[min(range(len(tau0)), key=lambda i: abs(tau0[i]-1))],'h(Tb=T')
print('1000 nm: Mean height of formation: ', hmean1, 'h(tau=1): ', h[min(range(len(tau1)), key=lambda i: abs(tau1[i]-1))])
print('1600 nm: Mean height of formation: ', hmean2, 'h(tau=1): ', h[min(range(len(tau2)), key=lambda i: abs(tau2[i]-1))])
print('5000 nm: Mean height of formation: ', hmean3, 'h(tau=1): ', h[min(range(len(tau3)), key=lambda i: abs(tau3[i]-1))])
Tb=func.BrightnessTemperature(5*1e-4,intt1*1e-4)
print(Tb)
#plt.figure(8)
#plt.plot(h,Tb)
#Tb1=func.BrightnessTemperature(1,contfunc1)
#Tb2=func.BrightnessTemperature(1.6,contfunc2)
#Tb3=func.BrightnessTemperature(5,contfunc3)
#print(Tb0)
#plt.figure(8)
#plt.plot(h,intt0)
#plt.plot(h,I)
#print(Tb0)
#print(temp)
#	plt.figure()
#	plt.title("H$^{-}$ transmittance at h = 0")
#	plt.plot(wavelength*1e-4,1/EXmin, color = "royalblue")
#	plt.grid(linestyle = "--")
#	plt.xlabel(r"Wavelength $\lambda$[$\mu$m]")
#	plt.ylabel(r"H$^{-}$ extinction$^{-1}$ [cm$^2$ per H-atom]$^{-1}$")
#	plt.subplots_adjust(bottom = 0.12)
##	plt.savefig(savepath1 + "extinction_transmittance.pdf")
#	plt.show()


#fig = plt.figure(4)
##temp = T(h=0), h[69]=0
#wav = wav*1e4
#nneutH=nhyd-nprot
#sigmaT= 6.648e-25
#EXmin = exthmin(0.5*1e4,temp,nel)*nneutH 
#Thomson = nel*sigmaT
#plt.semilogy(h,EXmin)
#plt.semilogy(h,Thomson)
#plt.semilogy(h,EXmin+Thomson)