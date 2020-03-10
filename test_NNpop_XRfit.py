import numpy as np
from NN_XRfit import *
import pylab as pl
import sys
sys.path.append('/home/mrinal/Backup/Python/Projects/Fortran_Routines/')
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
	
	
qmax=0.4
dir='/home/mrinal/Backup/UIC/Experiments/APS_Runs/2011dec/Extracted_Data/After_Run/'
##For reading experimental data file#
##----------------------------------#
fname='sample2_ref_0.45v.ref'
dat=np.loadtxt(dir+fname)
qoff=0.0
data=dat[np.argwhere(dat[:,0]<qmax)[:,0],:]
data[:,0]=data[:,0]-qoff
data[:,1]=data[:,1]#/np.average(data[0:1,1])
#data=np.delete(data,2,1)
##----------------------------------#

##For generating reflectivity curves#
##----------------------------------#
#q=np.arange(0,0.6,0.01)
#d=[0.0,30,40.0,10.0]
#rho=[0.0,0.6,0.25,0.334]
#beta=np.zeros_like(rho)
#sig=[0.0,3.0,3.0,3.0]
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
N=[30,30]
d=N#
rho=[0.7,0.38]
beta=np.zeros_like(rho)
rho_range=[0.16,1.0]
lamda=0.41328
rho_b=[0.334,0.38]
delta=rho_range[1]-rho_range[0]

maxiter=10**6
#data=np.vstack((data[:,0],data[:,1])).transpose()
mif=NN_XR(data,rho_b=rho_b,d=d,Nl=N,rho=rho,beta=beta,rho_range=rho_range,lamda=lamda,delta=delta/10,logy=0,Rf=1,Npop=10)
conv=int((2*np.pi/data[-1,0])*sum(N)/sum(d))
if np.mod(conv,2)==0:
  conv=conv+1
print conv
mif.evolve(qmax=qmax/3,iplot=1,iprint=1,maxiter=10000)
rho=[mif.rhol[i] for i in range(1,len(rho)+1)]
mif.create_layers(N,d,rho,beta)
mif.delta=delta/10
mif.create_pop(qmax=qmax)
#print [item['merit'] for item in mif.pop.values()]
#print mif.best_rho
mif.pop_evolve(qmax=qmax,iplot=10,iprint=10,deltam=0.001,avg=5)
#print mif.pop
#x=mif.x
#y=mif.y
#yerr=mif.y_err
##z=mif.z
q=np.arange(data[0,0],data[-1,0],0.001)
mif.rhol=mif.best_rho
fit=mif.func(q)
rho=mif.rhol
z=mif.z
rrf_fname=fname.split('.')[0]+'_rrf.txt'
fit_fname=fname.split('.')[0]+'_fit.txt'
rho_fname=fname.split('.')[0]+'_edp.txt'
np.savetxt(dir+rrf_fname,np.vstack((mif.x,mif.y,mif.y_err)).transpose())
np.savetxt(dir+fit_fname,np.vstack((q,fit)).transpose())
np.savetxt(dir+rho_fname,np.vstack((z,rho)).transpose())