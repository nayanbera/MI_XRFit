import numpy as np
import sys
import os
home=os.getenv('HOME')
sys.path.append(home+'/Documents/Python/Projects/MI_XRFit/')
from xr_ref import *
import pylab as pl
from scipy.special import gamma
from scipy.interpolate import UnivariateSpline, interp1d
import copy


class XRLSS:
    """
    Routine for Model Independent fitting of X-ray specular reflectivity
    """
    def __init__(self, data, datatype='xrr', rho_b=[0.0,0.334],beta_b=[0.0,0.0], d=[100.0],Nl=[50],rho=[0.5],rho_range=[0,1.0],beta=[0.0],deltarho=1.0,lamda=1.54,T=22.0, gam=73.0,lsd=560,s1=0.020,thetain=0.0904,qpar=0.01,qmax=1.3,bkg=0.0,logy=0,Npop=1,iplot=1):
        """
        Class to do model independent fitting of X-ray reflectivity data.
        data=2 column oar 3 column reflectivity data with the sequence of columns [q, ref[, ref_err]]
        datatype='xrr' for specular reflectivity or 'gix' gixos data 
        rho_b=Electron densities in el/Angstorms^3 of the bulk phases
        beta_b=absorption coefficients of the bulk phases
        d=possible film thickness of the layers in Angstroms
        rho=list of possible electron densities of each layer.
        N=list of Number of sublayers in each layer of the film
        rho_range=possible range of electron densities of the film in the form of [min,max]
        beta=list of absorption coefficients of the layers
        deltarho=maximum change allowed in the electron density of the layers.
        lamda=wavelengths of x-ray used
        T=Temperature of the third phase in degree celcius
        gam=interfacial tension of the air-water interface in mN/m
        qpar=in-plave wavevector in Ang^-1 at which GIXOS is measured
        lsd=sample to detector distance in mm
        s1=vertical incident slit in mm
        thetain=inicident angle in degrees
        qmax=Maximum in-pane wave-vector corresponding to the system i.e 2pi/d for water with d being the diameter of water molecule
        logy=0 means data and 1 means log(data) will be used for fitting
        Rf= 0 means data and 1 means data/rf will be used
        iplot=1 means plot or 0 means do not plot
        """
        self.lamda=lamda
        self.Nlayer=np.array(Nl)
        self.rho_b=np.array(rho_b)
        self.beta_b=np.array(beta_b)
        self.d=np.array(d)
        self.rho=np.array(rho)
        self.beta=np.array(beta)
        self.T=T
        self.gam=gam    
        self.qpar=qpar    
        self.rho_range=rho_range
        self.Npop=Npop
        self.data=data
        self.datatype=datatype
        self.create_layers(self.Nlayer,self.d,self.rho,self.beta)
        self.qc=4.0*np.pi*np.sin(np.sqrt(2.817e-5*np.average(self.rhol)/np.pi))
        self.norm=0#np.argwhere(self.data[:,0]<self.qc)[-1][0]
        self.prepare_data()
        self.deltarho=deltarho
        self.logy=logy
        self.qmax=qmax
        self.bkg=bkg
        #self.qzmax,self.qzmin=self.footprint(self.x,lsd,s1,thetain)
        #N=np.ones_like(self.Nlayer)
        self.z=np.cumsum(self.dl)
        #self.z[-1]=self.z[-2]+self.dl[-2]
        if iplot!=0:
            self.fig=pl.figure(figsize=(5,9))            
            self.ax1=self.fig.add_subplot(311)
            if self.logy!=0:
                self.ax1.set_yscale('log')
            self.ax1.errorbar(self.x,self.y,self.y_err,fmt='ko')
            self.ax3=self.fig.add_subplot(312)
            self.ax3.plot(self.x,(self.y-self.func(self.x))/self.y_err,'g-')
            self.chimax=np.max(np.abs((self.y-self.func(self.x))/self.y_err))
            self.ax3.set_ylim(-self.chimax,self.chimax)
            #self.q=np.arange(self.x[0],self.x[-1],0.001)
            self.ax1.plot(self.x,self.func(self.x),'r-')
            self.ax2=self.fig.add_subplot(313)
            #width=np.append(self.dl[1:],[0.0],0)
            #print self.z,width
            #self.ax2.bar(self.z,(np.array(self.rhol)-self.rho_b[0])/(self.rho_b[1]-self.rho_b[0]),width=self.dl)
            self.ax2.step(self.z,(np.array(self.rhol)-self.rho_b[0]),'b',lw=2)#/(self.rho_b[1]-self.rho_b[0]),'b-',)
            self.ax3.set_xlabel(r'Q$_z$($\AA^{-1}$)')
            self.ax3.set_ylabel(r'(I-I$_{fit}$)/I$_{err}$')
            self.ax1.set_xlabel(r'Q$_z$($\AA^{-1}$)')
            self.ax1.set_ylabel(r'I')
            self.ax2.set_xlabel(r'z ($\AA$)')
            self.ax2.set_ylabel(r'$\rho$ (el/$\AA^{3}$)')
            pl.tight_layout()
            #self.fig.canvas.draw()
            self.fig.show()
    
    def prepare_data(self):
        """
        Prepare the data for the fitting
        """
        self.data=np.array([dat for dat in self.data if dat[-1]>0])
        try:
            self.x=self.data[:,0]
            if self.datatype=='xrr':
                rf,r=parratt(self.x,self.lamda,[0.0,0.0],self.rho_b,self.beta_b)
                try:
                    self.y=self.data[:,1]/rf
                    self.y_err=self.y*0.005#self.data[:,2]/rf                    
                except:
                    print "Errorbars are created with 2% of data!!!"
                    self.y_err=0.005*self.y #0.5% errorbar w.r.t that data if errorbars are not provided
    
            else:
                try:
                    self.y=self.data[:,1]*self.data[:,0]**2/self.data[self.norm,1]/self.data[self.norm,0]**2 # Normalizing the data with the first point for gixos
                    self.y_err=self.data[:,2]*self.data[:,0]**2/self.data[self.norm,1]/self.data[self.norm,0]**2
                except:
                    print "Errorbars are created with 2% of data!!!"
                    self.y_err=0.02*self.data[:,1]*self.data[:,0]**2/self.data[self.norm,1]/self.data[self.norm,0]**2 #2% errorbar w.r.t that data if errorbars are not provided            
        except:
            print "Error:: Please provide valid data!!"
    
    
    def create_layers(self,N,d,rho,beta):
        """
        Create the layers for the fitting
        N=list of number of sublayers in each layer.
        d=list of layer thicknesses.
        rho=list of the electron densities of the layers
        beta=list of the absorption coefficients of the layers
        """
        self.rhol=[self.rho_b[0]]
        self.betal=[self.beta_b[0]]
        self.dl=[0.0]
        for j in range(len(N)):
            for i in range(1,N[j]+1):
                self.rhol.append(rho[j])
                self.betal.append(beta[j])
                self.dl.append(d[j]/N[j])
        self.rhol.append(self.rho_b[1])
        self.betal.append(self.beta_b[1])
        self.dl.append(self.d[0]/2)#(d[-1]/N[-1])
        self.rhol=np.array(self.rhol)
        self.betal=np.array(self.betal)
        self.dl=np.array(self.dl)

    def create_pop(self,qmax=0.1):
        """
        Create the population of Electron Density profiles
        """
        self.pop={}
        for i in range(self.Npop-1):
            np.random.seed()
            #rho=np.random.rand(len(self.rho))*(self.rho_range[1]-self.rho_range[0])
            self.create_layers(self.Nlayer,self.d,self.rho,self.beta)
            self.pop[i]={'rho':self.rhol,'merit':self.merit(qmax)}
        self.create_layers(self.Nlayer,self.d,self.rho,self.beta)
        self.pop[self.Npop-1]={'rho':self.rhol,'merit':self.merit(qmax)}
        self.pop_stats()
    
        
    def smooth_pop(self):
    #self.delta=self.delta*(1.0-1e-3)
        for key in self.pop.keys():
            tmprho=self.pop[key]['rho']
            tnlayer=len(tmprho)-2
            for k in range(tnlayer,0,-1):
                srho=0.0
                nr=0
                for kc in range(-1,2,2):
                    if tnlayer+2>k+kc>-1:
                        srho=srho+tmprho[k+kc]
                        nr=nr+1
                    elif k+kc<=-1:
                        srho=srho+self.rho_b[0]
                        nr=nr+1
                    elif tnlayer+2<=k+kc:
                        srho=srho+self.rho_b[1]
                        nr=nr+1
                self.pop[key]['rho'][k]=srho/nr
        self.pop_stats()
    
    def minmax_rho(self,rho,i):
        """
        Calculates maximum and minimum possible value of rho of ith layer keeping in mind the rho of (i+1)th and (i-1)th layer
        rho=list of electron densities of all the sublayers.
        i=ith layer
        delta=maximum change in the electron density of ith layer
        """
        np.random.seed()
        if rho[i+1]!=rho[i-1]:
            return np.min([np.max([rho[i+1],rho[i-1]])+self.deltarho*np.random.rand(),self.rho_range[1]]), np.max([np.min([rho[i+1],rho[i-1]])-self.deltarho*np.random.rand(),self.rho_range[0]])
            #return np.min([rho[i+1],rho[i-1]]),np.max([rho[i+1],rho[i-1]])
        else:
            return np.min([rho[i]+self.deltarho*np.random.rand(),self.rho_range[1]]),np.max([self.rho_range[0],rho[i]-self.deltarho*np.random.rand()])
    
        
    def pop_evolve(self,qmax=0.1,deltam=1,maxiter=10**6,iplot=0,iprint=0,avg=0):
        """
        Evolve the population
        """
        niter=0
        bestm=self.best_merit
        while niter<=maxiter and self.delta>deltam:
            try:
                for i in self.pop.keys():
                    self.rhol=copy.copy(self.pop[i]['rho'])
                    self.evolve(qmax=qmax,maxiter=1,deltam=deltam)
                    self.pop[i]['rho']=copy.copy(self.rhol)
                    self.pop[i]['merit']=self.merit(qmax)
                niter=niter+1
                self.pop_stats()
                self.delta=self.delta*(1.0-1e-4)
                if avg!=0 and np.mod(niter,avg)==0 and self.best_merit<bestm:
                    bestm=self.best_merit
                    self.smooth_pop()
                if iprint!=0 and np.mod(niter,iprint)==0:
                    self.rhol=self.mean_rho
                    merit=self.merit(qmax)
                    print niter, 'merit= ', bestm, 'delta= ',self.delta#,'rho= ', self.mean_rho,'+/- ', self.std_rho
                if iplot!=0 and np.mod(niter,iplot)==0:
                    self.plot_pop()            
            except KeyboardInterrupt:
                break
        
    def pop_stats(self):
        """
        Calculates statistics about the population
        """
        self.best_merit=self.pop[0]['merit']
        self.best_rho=self.pop[0]['rho']
        rho=np.zeros_like(self.rhol)
        for pop in self.pop.values():
            if pop['merit']<self.best_merit:
                self.best_merit=pop['merit']
                self.best_rho=pop['rho']
            rho=np.vstack((rho,pop['rho']))
        rho=np.delete(rho,0,axis=0)
        self.mean_rho=np.mean(rho,axis=0)
        self.std_rho=np.std(rho,axis=0)
        self.min_rho=np.amin(rho,axis=0)
        self.max_rho=np.amax(rho,axis=0)

     
    def plot_pop(self):
        """
        Plotting the results
        """
        self.ax1.clear()
        self.ax2.clear()
        self.ax3.clear()
        self.ax1.errorbar(self.x,self.y,self.y_err,fmt='b.')
        merits=[pop['merit'] for pop in self.pop.values()]
        self.z=np.cumsum(self.dl)
        self.ax3.semilogy(merits,'r*')
        #for pop in self.pop.values():
            #self.rhol=pop['rho']
            #self.ax1.plot(self.q,self.func(self.q),'r-')
            #self.ax2.plot(self.z,(np.array(self.rhol)-self.rho_b[0])/(self.rho_b[1]-self.rho_b[0]),'r-')
            #self.ax2.plot(self.z,self.rhol,'r-')
        self.rhol=self.mean_rho#self.best_rho
        self.ax1.plot(self.q,self.func(self.q),'b-')
        self.rhol=copy.copy(self.best_rho)
        self.ax1.plot(self.q,self.func(self.q),'r-')
        #self.ax2.plot(self.z,(np.array(self.rhol)-self.rho_b[0])/(self.rho_b[1]-self.rho_b[0]),'b-')
        #self.ax2.plot(self.z,self.best_rho,'r-')
        #self.ax2.errorbar(self.z,self.mean_rho,self.std_rho,fmt='b-')
        self.ax2.bar(self.z,(np.array(self.rhol)-self.rho_b[0])/(self.rho_b[1]-self.rho_b[0]),width=self.dl)
        #self.rhol=self.mean_rho
        #self.ax1.plot(self.q,self.func(self.q),'g-')
        #self.ax2.plot(self.z,(np.array(self.rhol)-self.rho_b[0])/(self.rho_b[1]-self.rho_b[0]),'g-')
        self.fig.canvas.draw()
    
    
    
    def evolve(self,qmax=0.1,deltam=1,maxiter=10**6,iplot=0,iprint=0,avg=0,Navg=0):
        """
        Evolve the electron densities of the layers 
        qmax = maximum q value of the data to be used for the fitting
        deltam=minimum change in merit function or chi square value for stopping the iterations
        iplot=0 or N for not or plotting the results after N iterations
        iprint=0 or N for not or printing the results after N iterations
        avg=Number of sublayers to be used for averaging
        Navg=After how many iterations averging will be done. Navg=0 means no avering will be done 
         """
        if avg>0:
            win=np.ones(avg)/avg
        bestm=self.merit(qmax)
        bestrho=copy.copy(self.rhol)
        bestd=copy.copy(self.dl)
        bestbeta=copy.copy(self.betal)
        dm=bestm
        iterations=0
        tnlayer=len(self.rhol)-2
        chisq=np.array([0,bestm])
        bestbestm=bestm
        bestrho=copy.copy(self.rhol)
        change=0
        nochange=0
            
        while iterations<maxiter and self.deltarho>0.001:
            try:  
                for j in range(1,tnlayer+1):
                    if tnlayer>1:
                        k=np.random.random_integers(1,tnlayer)
                    else:
                        k=1
                    maxr,minr=self.minmax_rho(self.rhol, k)
                    if abs(maxr-minr)>0.01:
                        np.random.seed()
                        temp_rho=minr+np.random.rand()*(maxr-minr)
                        if self.rho_range[1]>temp_rho>self.rho_range[0]:
                            old_rho=copy.copy(self.rhol[k])
                            self.rhol[k]=temp_rho
                            m1=self.merit(qmax)
                            if m1>=bestm:
                                self.rhol[k]=copy.copy(old_rho)
                            else:
                                bestm=m1
                                change=change+1
                        while np.abs(self.rhol[k]-self.rho_range[0])<deltam:
                            self.rhol=np.delete(self.rhol,1,0)
                            self.rhol=np.append(self.rhol,[self.rho_b[1]],0)
                        bestm=self.merit(qmax)

                if bestm<bestbestm:
                    if np.abs(bestm-bestbestm)<0.00001:
                        nochange=nochange+1
                        bestbestm=copy.copy(bestm)
                        bestrho=copy.copy(self.rhol)
                        pl=1
                        if nochange>10:
                            break
                    else:
                        nochange=0
                    bestbestm=copy.copy(bestm)
                    bestrho=copy.copy(self.rhol)
                    pl=1
                else:
                    pl=0
                tmprho=copy.copy(self.rhol)
                sumrho=np.sum(tmprho[1:-1])
                if Navg!=0 and np.mod(iterations,Navg)==0:
                    self.rhol[1:-1]=np.convolve(self.rhol,win,mode='valid')
                    change=0
                    if self.merit(qmax)<2.0*bestm:
                        self.deltarho=self.deltarho*(1.0-1e-1)                        
                    else:
                        self.rhol=copy.copy(bestrho)
                bestm=self.merit(qmax)
                iterations=iterations+1
                if iprint!=0:# and pl!=0:
                    print iterations,bestm,bestbestm,self.deltarho,tnlayer,nochange
                trho=copy.copy(self.rhol)
                if iplot!=0 and pl!=0:
                    self.ax1.clear()
                    self.ax2.clear()
                    self.ax3.clear()
                    if self.logy!=0:
                        self.ax1.set_yscale('log')
                    chisq=np.vstack((chisq,np.array([iterations,bestbestm])))
                    self.z=np.cumsum(self.dl)
                    width=np.append(self.dl[1:],[0.0],0)
                    #self.z[-1]=self.z[-2]+self.dl[-2]
                    self.ax1.errorbar(self.x,self.y,self.y_err,fmt='r.',ecolor='red')
                    self.ax1.axvline(x=qmax,color='red')
                    self.ax1.set_xlim(np.min(self.x),np.max(self.x))
                    #print chisq
                    self.ax3.semilogy(chisq[:,0],chisq[:,1],'r*')
                    q=np.arange(self.x[0],self.x[-1],0.001)
                    self.ax1.plot(q,self.func(q),'b-',lw=2)
                    self.ax2.step(self.z,self.rhol,'r-')
                    self.rhol=bestrho
                    self.ax1.plot(q,self.func(q),'g-',lw=2)
                    self.ax2.step(self.z,self.rhol,'g-')
                    self.fig.tight_layout()
                    self.fig.canvas.draw()
                    self.fig.canvas.flush_events()
                    if np.abs(chisq[-1,1]-np.average(chisq[-10:-1,1]))*100/chisq[-1,1]<1.0:
                        break
                self.rhol=trho
            except (KeyboardInterrupt, SystemExit):
                self.rhol=copy.copy(bestrho)
                break
        print 'Chisq= ',self.merit(qmax)
        print "Fitting completed!!"


    def cal_t(self,qz):
        #tsc=np.abs(2*qz/(qz+np.sqrt(qz**2-qc**2+0j)))**2
            ref,trans=parratt_mat(qz,self.lamda,self.dl,self.rhol,self.betal)
            tsc=np.abs(trans)**2
            #print parratt_mat(qz,self.lamda,self.dl,self.rhol,self.betal)
            return tsc
    
    def func(self,x):
        """
        Calculates the fitting function at the q values.
        x=list of q values
        """
        if self.datatype=='xrr':
            refq,r=parratt(x,self.lamda,self.dl,self.rhol,self.betal)
            rf,r=parratt(x,self.lamda,[0.0,0.0],self.rho_b,self.beta_b)
            fy=(refq+self.bkg)/rf
        else:
            xnorm=x[0]#np.argwhere(x<self.qc)[-1]
            fy=self.gixos(x)
            fy=fy/fy[xnorm]+self.bkg
        return fy
    
    def phiqz(self,x):
        """
        Calculates structure factor for GIXOS calculation
        """
        if type(x)!=list:
            x=[x]
        z=np.cumsum(self.dl)
        drho=np.diff(self.rhol)
        phi=[]
        for q in x:
            phi.append(np.sum(drho*np.exp((0+1j)*q*z[:-1])))
        return np.abs(np.array(phi))**2

    def capillary(self,qz,qpar,gam,T,qmax):
        """
        Calculates the diffuse scattering part of scattering due to capillary wave
        """
        KbT=1.38e-23*(273.0+T)
        eta=KbT*qz**2*1e23/2.0/np.pi/gam
        return 2**(1-eta)*gamma(1.0-eta/2.0)/gamma(eta/2.0)/qpar**(2-eta)/qz**3
        #return (qpar*np.exp(-0.5772)/qmax)**eta*eta/qpar**2/qz**2

    def gixos(self,qz):
        """
        Calculates the gixos as a funciton of qz at a particular qpar
        """
        gix=[]
        #self.qc=4.0*np.pi*np.sin(np.sqrt(2.817e-5*np.average(self.rhol)/np.pi))
        for qz1 in qz:
            #fsum=0.0
            #for qz1 in [qz[i]]:#np.linspace(self.qzmin[i],self.qzmax[i],10):
            fsum=self.capillary(qz1,self.qpar,self.gam,self.T,self.qmax)*self.phiqz([qz1])[0]
            gix.append(fsum)#*self.cal_t(qz1)[0])
        return np.array(gix)*qz**2

    def footprint(self,qz,lsd,s1,thetain):
        """
        Apply the resolution function on the calculated scattered intensity due to extended footprint 
        of the x-rays where
        qz=reciprocal wave-vector along z direction in Angs^-1
        lsd=sample to detector distance in mm
        s1=vertical incident slit in mm
        thetain=incident angle in degrees
        """
        lamda=self.lamda
        theta=np.arcsin(lamda*qz/2.0/np.pi-np.sin(thetain*np.pi/180.0))
        l1=s1/np.sin(thetain*np.pi/180.0)
        return 2*np.pi*(np.sin(np.arctan(lsd*np.tan(theta)/(lsd-l1/2)))+np.sin(thetain*np.pi/180.0))/lamda\
        ,2*np.pi*(np.sin(np.arctan(lsd*np.tan(theta)/(lsd+l1/2)))+np.sin(thetain*np.pi/180.0))/lamda

      
    def merit(self,qmax):
        """
        Calculates the merit function or Chi square upto the qvalue given by qmax
        """
        xl=self.x[np.argwhere(self.x<=qmax)[:,0]]
        xil=np.sum((self.y[:len(xl)]-self.func(xl))**2/self.y_err[:len(xl)]**2)/len(xl)#*np.sum(np.sqrt(1.0+np.diff(self.rhol)**2))
        #if len(xr)>0:
        #    xir=np.sum((self.y[:len(xr)]-self.func(xr))**2/self.y_err[:len(xr)]**2)
        #else:
        #    xir=0.0
        return xil#+0.1*xir

    
        
        
