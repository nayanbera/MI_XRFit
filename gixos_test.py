import numpy as np
import sys
import os
import pylab as pl
home=os.getenv('HOME')
sys.path.append(home+'/Documents/Python/Projects/MI_XRFit')
data_dir=home+'/Documents/Binhua_paper/XR_raw/'
from MI_XRLSS import XRLSS
data=np.loadtxt(data_dir+'ref_295-304_raw.txt')
ds=[140]#,40]
rho=[0.334]#,0.7]
beta=np.zeros_like(rho)
N=[4]
lamda=1.24
start=0
mif=XRLSS(data[start:,:],datatype='xrr',rho_b=[0.0,0.334],d=ds,rho=rho,deltarho=0.5,beta=beta,Nl=N,lamda=lamda,T=22.0, gam=72.3-22.8,qpar=0.01,iplot=1,rho_range=[0.0,2.0],logy=1)
while mif.merit(0.1)>1000:
	mif.create_layers(N,ds,rho,beta)
	mif.evolve(qmax=0.1,deltam=0.1,maxiter=1000,iplot=1,iprint=1,avg=0)
	print mif.merit(0.1)
rho=mif.rhol[1:-1]
N=N*N[0]
d=mif.dl[1:-1]
beta=mif.betal[1:-1]
mif.deltarho=0.25
mif.create_layers(N,d,rho,beta)
mif.evolve(qmax=0.2,deltam=0.05,maxiter=1000,iplot=1,iprint=1,avg=3)
rho=mif.rhol[1:-1]
N=N*N[0]
d=mif.dl[1:-1]
beta=mif.betal[1:-1]
mif.deltarho=0.125
mif.create_layers(N,d,rho,beta)
mif.evolve(qmax=0.4,deltam=0.025,maxiter=1000,iplot=1,iprint=1,avg=3)
rho=mif.rhol[1:-1]
N=N*N[0]
d=mif.dl[1:-1]
beta=mif.betal[1:-1]
mif.deltarho=0.125
mif.create_layers(N,d,rho,beta)
mif.evolve(qmax=0.8,deltam=0.025,maxiter=1000,iplot=1,iprint=1,avg=3)
rhol=mif.rhol
z=mif.z
x=mif.x
y=mif.y
err=mif.y_err
fit=mif.func(x)
print len(mif.z),len(mif.rhol)
pl.savetxt(data_dir+'aed_349-356.txt',np.vstack((z,rhol)).transpose())
pl.savetxt(data_dir+'fit_349-356.txt',np.vstack((x,y,err,fit)).transpose())
#pl.show()
