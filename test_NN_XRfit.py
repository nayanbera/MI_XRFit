import numpy as np
from NN_XRfit import *
import pylab as pl
import sys
from xr_ref import *

def edp_gen(d=[0.0,0.0],rho=[0.0,0.334],sig=[0.0,3.0]):
    z=np.arange(-5*sig[1],5*sig[-1]+1.0+np.sum(d),1.0)
    drho=np.zeros_like(z)
    zl=0.0
    for i in range(1,len(d)):
        drho=drho+(rho[i]-rho[i-1])*np.exp(-(z-zl)**2/2/sig[i]**2)/2.507/sig[i]
        zl=zl+d[i]
    srho=np.cumsum(drho)
    return z,drho,srho


qmax=0.6
#dir='C:/Users/Mrinal Bera/Documents/Collaborations/Pan/Iron_Nanoparticles/'
dir='./refdata/MDX1342/'
###For reading experimental data file#
###----------------------------------#
#fname='FeO_NP_rrf.txt'
fname='s_63_s67_MDX1342_0.5mgml_rrf.txt'

#dir='./refdata/'
##For reading experimental data file#
##----------------------------------#
#fname='s159_s163_PS80_0.2CMC_rrf.txt'

data=np.loadtxt(dir+fname)
#data=np.delete(data,2,1)
##----------------------------------#

##For generating reflectivity curves#
##----------------------------------#
#q=np.arange(0,0.6,0.01)
#d=[0.0,20,40.0,10.0]
#rho=[0.0,0.6,0.25,0.334]
#beta=np.zeros_like(rho)
#sig=[0.0,5.0,5.0,3.0]
#dat,r=parratt_born(q,1.24,d,rho,beta,sig)
#z1,drho,srho=edp_gen(d=d,rho=rho,sig=sig)
#data=np.vstack((q,dat)).transpose()
#fname='ref_3.txt'
#np.savetxt(dir+fname,data)
#ax=pl.figure()
#ax.add_subplot(211)
#pl.semilogy(q,dat,'ro')
#ax.add_subplot(212)
#pl.plot(z1,srho,'r-')
#pl.show()
#gedp=np.vstack((z1,srho)).transpose()
#np.savetxt(dir+fname.split('.')[0]+'_gedp.txt',gedp)
##-----------------------------------#

##For giving the parameters for fitting#
##-------------------------------------#
"Other guess values- By Mrinal"
#N=np.ones(10)
#d=np.ones(10)*18#[20,137]#
#rho=np.ones(10)*0.34
#rho[-1]=0.42
#beta=np.zeros_like(rho)
#rho_range=[0,2.0]
#lamda=1.24
#rho_b=[0.0,0.334]
#delta=rho_range[1]-rho_range[0]


N=[10]
d=[180]
rho=[0.4]
#rho[-1]=0.42
beta=np.zeros_like(rho)
rho_range=[0,2.0]
lamda=1.24
rho_b=[0.0,0.334]
delta=rho_range[1]-rho_range[0]




maxiter=10**6
#data=np.vstack((data[:,0],data[:,1])).transpose()
mif=NN_XR(data,rho_b=rho_b,beta_b=[0.0,0.0],d=d,Nl=N,rho=rho,beta=beta,rho_range=rho_range,lamda=1.24,delta=delta,logy=0,Rf=1)

mif.evolve(qmax=qmax/3,iplot=20,deltam=0.1,iprint=3,maxiter=10000,saveFig=True,dir=dir+'Images/')
rho=[mif.rhol[i] for i in range(1,len(mif.rhol)-1)]
beta=[mif.betal[i] for i in range(1,len(mif.betal)-1)]
d=[mif.dl[i] for i in range(1,len(mif.dl)-1)]
print(len(rho),len(beta),len(d))
N=np.ones_like(rho)*10
mif.create_layers(N,d,rho,beta)
N=len(rho)
print (N)
#mif.delta=0.05
#avg : eaverages 3 slabs and gives a value, always an odd number. Average of 2 slabs is not acceptable
mif.evolve(qmax=qmax,iplot=10,iprint=3,deltam=0.01,avg=3,saveFig=True,dir=dir+'Images/')
rho=[mif.rhol[i] for i in range(1,len(mif.rhol)-1)]
beta=[mif.betal[i] for i in range(1,len(mif.betal)-1)]
d=[mif.dl[i] for i in range(1,len(mif.dl)-1)]
print(len(rho),len(beta),len(d))
N=np.ones_like(rho)*2
mif.create_layers(N,d,rho,beta)
N=len(rho)
print(N)
mif.evolve(qmax=qmax,iplot=10,iprint=3,deltam=0.001,avg=3,saveFig=True,dir=dir+'Images/')
x=mif.x
y=mif.y
yerr=mif.y_err
#z=mif.z
q=np.arange(data[0,0],data[-1,0],0.001)
fit=mif.func(q)
rho=mif.rhol
z=mif.z
print(len(rho),len(z))
#fit_fname=fname.split('.')[0]+'_fit1.txt'
fit_fname=fname+'_fit1.txt'
#rho_fname=fname.split('.')[0]+'_edp1.txt'
rho_fname=fname+'_edp1.txt'
edp=np.vstack((z,rho)).transpose()
edp=pl.insert(edp,0,[z[0]-20,rho[0]],axis=0)
edp[-1,0]=z[-1]+20
np.savetxt(dir+fit_fname,np.vstack((q,fit)).transpose())
np.savetxt(dir+rho_fname,edp)
