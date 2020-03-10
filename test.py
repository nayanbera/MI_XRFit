import numpy as np
import pylab as pl
import xr_ref as xr

q=np.arange(0,0.7,0.001) #qvalues
lam=1.24 #Wavelength of x-rays used
d=[0.0,100.0,1e4] #List of Thicknesses of the layers. Here only 3 layers are used.
rho=[0.0,0.42,0.71] #List of Electron densities of the layers. 
beta=[0.0,3.e-7,4e-7] #List of beta values of the layers where refractive index is n=1-delta-i*beta. Here no absorption is taken.
sig=[0.0,4.0,4.0] #List of roughnesses

ref1,ref2=xr.parratt_born(q,lam,d,rho,beta,sig)
pl.semilogy(q,ref1,'g-',lw=2)
pl.show()