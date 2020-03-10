import numpy as np
import sys
sys.path.append('/home/mrinal/Backup/Python/Projects/Fortran_Routines/')
from xr_ref import *
import pylab as pl
import copy

class NN_XR:
    """
    Routine for Model Independent fitting of X-ray specular reflectivity
    """
    def __init__(self, data, rho_b=[0.0,0.334],beta_b=[0.0,0.0], d=[100.0],Nl=[50],rho=[0.5],rho_range=[0,1.0],beta=[0.0],delta=1.0,lamda=1.54,T=22.0, gam=73.0,logy=0,Rf=0,Npop=1,iplot=1):
        """
        Class to do model independent fitting of X-ray reflectivity data.
        Input: 
        data=2 column or 3 column reflectivity data with the sequence of columns [q, ref[, ref_err]]
        rho_b=Electron densities in el/Angstorms^3 of the bulk phases
        beta_b=absorption coefficients of the bulk phases
        d=possible film thickness of the layers in Angstroms
        rho=list of possible electron densities of each layer.
        N=list of Number of sublayers in each layer of the film
        rho_range=possible range of electron densities of the film in the form of [min,max]
        beta=list of absorption coefficients of the layers
        delta=
        deltam = points upto which the chi square value should minimize to , say deltam = 0.001
        lamda=wavelengths of x-ray used
        T=Temperature of the third phase in degree celcius
        gam=interfacial tension of the air-water interface in mN/m
        logy=0 means data and 1 means log(data) will be used for fitting
        Rf= 0 means data and 1 means data/rf will be used
        iplot=1 means plot or 0 means do not plot
        """
        
        """This value can change depending on when you call mif.evolve. So you start with say N = 10 fo qmax=qmax/3
        Then you call mif.evolve for the same (refer test_ file). Then you divide that each of the N layers 
        to further N sublayers and call mif.evolve again for the whole range of q"""
        self.lamda=lamda
        self.Nlayer=Nl                                                         
        self.rho_b=rho_b
        self.beta_b=beta_b
        self.d=d
        self.rho=rho
        self.beta=beta
        
        self.rho_range=rho_range
        self.Npop=Npop
        self.prepare_data(data=data,logy=logy,Rf=Rf)
        self.delta=delta
        self.Rf=Rf
        self.logy=logy
        #N=np.ones_like(self.Nlayer)
        self.create_layers(self.Nlayer,self.d,self.rho,self.beta)
        self.z=np.cumsum(self.dl)
        self.counter=1
        #self.z[-1]=self.z[-2]+self.dl[-2]
        if iplot!=0:
            self.fig=pl.figure(figsize=(8,6))
            self.gs=self.fig.add_gridspec(2,2)
            self.ax1=self.fig.add_subplot(self.gs[0,0])
            self.ax1.set_title('Fresnel Normalized Reflectivity')
            self.ax1.set_yscale('log')
            self.ax1.errorbar(self.x,self.y,self.y_err,fmt='b.')               #Plot is for R/Rf vs qz with error, x and y comes from class prepare_data"""

            self.ax3=self.fig.add_subplot(self.gs[1,:])
            self.ax3.set_title('Fitting Info')
            #self.ax3.plot(self.x,(self.y-self.func(self.x))/self.y_err,'g-')
            self.q=np.arange(self.x[0],self.x[-1],0.001)                       #Range of qz"""
            self.ax1.plot(self.q,self.func(self.q),'r-')                       #This plot is in RED , it tries to fin the minima other than the GREEN curve, if it finds it it accepts
                                                                               #it else goes to GREEN curve"""
            self.ax2=self.fig.add_subplot(self.gs[0,1])                                 #Plots the edp"""
            self.ax2.set_title('Electron Density Profile')
            self.ax2.set_xlabel('Depth (\212B)')
            self.ax2.set_ylabel(r'$\Chi^2$')
            #width=np.append(self.dl[1:],[0.0],0)
            #print self.z,width
            #self.ax2.bar(self.z,(np.array(self.rhol)-self.rho_b[0])/(self.rho_b[1]-self.rho_b[0]),width=self.dl)
            #pl.show()
    
    def prepare_data(self,data=[],logy=0,Rf=0):
        """
        Prepare the data for the fitting
        x = qz and y = R/Rf
        logy=0 or 1 for linear or logarthmic of R/Rf before using it for Chi square calculaiton
        Rf=0 or 1 for dividing the experimental data with fresnel reflectivity before using it for Chi square calculation
        """
        try:
            self.x=data[:,0]
            self.y=data[:,1]
            rf,r=parratt(self.x,self.lamda,[0.0,0.0],self.rho_b,self.beta_b)
            try:
                self.y_err=data[:,2]
            except:
                self.y_err=1e-3*np.ones(len(self.x))
                if logy!=0:
                    """Converts into log plot using 2.303"""
                    self.y_err=self.y_err/self.y/2.303
                    self.y=np.log10(self.y)    
        except:
            print ("Error:: Please provide valid data!!")
    
    
    def create_layers(self,N,d,rho,beta):
        """
        Create the layers for the fitting, Nl, betal, rhol
        N=list of number of sublayers in each layer.
        d=list of layer thicknesses.
        rho=list of the electron densities of the layers
        beta=list of the absorption coefficients of the layers
        """
        self.rhol=[self.rho_b[0]]
        self.betal=[self.beta_b[0]]
        self.dl=[0.0]
        for j in range(len(N)):
            for i in np.arange(1,N[j]+1):
                self.rhol.append(rho[j])
                self.betal.append(beta[j])
                self.dl.append(d[j]/N[j])
        self.rhol.append(self.rho_b[1])
        self.betal.append(self.beta_b[1])
        self.dl.append(d[-1]/N[-1])

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
    
    def minmax_rho(self,rho,i,delta):
        """
        Calculates maximum and minimum possible value of rho of ith layer keeping in mind the rho of (i+1)th and (i-1)th layer
        rho=list of electron densities of all the sublayers.
        i=ith layer
        delta=maximum change in the electron density of ith layer
        """
        np.random.seed()
        if rho[i+1]!=rho[i-1]:
            #return rho[i-1],rho[i+1]
            return np.max([rho[i+1],rho[i-1]])+np.random.rand()*delta, np.min([rho[i+1],rho[i-1]])-delta*np.random.rand()
        else:
            #return rho[i+1],rho[i-1]
            return rho[i]+delta,rho[i]-delta*np.random.rand()
        #if rho[i+1]>rho[i-1]:
        #    return rho[i-1]+delta, rho[i+1]-delta
        #else:
        #    return rho[i+1]+delta, rho[i-1]-delta 
    
        
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
                    print (niter, 'merit= ', bestm, 'delta= ',self.delta)#,'rho= ', self.mean_rho,'+/- ', self.std_rho)
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
        self.ax2.plot(self.z,self.best_rho,'r-')
        self.ax2.errorbar(self.z,self.mean_rho,self.std_rho,fmt='b-')
        #self.rhol=self.mean_rho
        #self.ax1.plot(self.q,self.func(self.q),'g-')
        #self.ax2.plot(self.z,(np.array(self.rhol)-self.rho_b[0])/(self.rho_b[1]-self.rho_b[0]),'g-')
        self.fig.canvas.draw()
    
    
    
    def evolve(self,qmax=0.1,deltam=1,maxiter=10**6,iplot=0,iprint=0,avg=0,saveFig=False,dir=None):
        """
        Evolve the electron densities of the layers 
        qmax =maximum q value of the data to be used for the fitting
        deltam=the maximum change allowed in the electron density of an individual layer
        iplot=0 or N for not or plotting the results after N iterations
        iprint=0 or N for not or printing the results after N iterations
        avg=Number of sublayers to be used for averaging
        saveFig = True or False for saving the plots as images to make video later on
        dir = directory to save the figures
         """
        bestm=self.merit(qmax)
        bestrho=copy.copy(self.rhol)
        bestd=copy.copy(self.dl)
        bestbeta=copy.copy(self.betal)
        dm=bestm
        iterations=0
        tnlayer=len(self.rhol)-2
        chisq=np.array([1,bestm])
        bestbestm=bestm
        
        try:
            repeat=0
            while iterations<maxiter and self.delta>deltam:
                for j in range(1,tnlayer+1):
                    if tnlayer>1:
                        k=np.random.randint(1,tnlayer)
                    else:
                        k=1
                    maxr,minr=self.minmax_rho(self.rhol, k, self.delta)
                    if abs(maxr-minr)>1e-5:
                        np.random.seed()
                        temp_rho=minr+np.random.rand()*(maxr-minr)
                        if self.rho_range[1]>temp_rho>self.rho_range[0]:
                            old_rho=self.rhol[k]
                            self.rhol[k]=temp_rho
                            m1=self.merit(qmax)
                            if m1>=bestm:
                                self.rhol[k]=old_rho
                            else:
                                bestm=m1
                    while np.abs(self.rhol[1]-self.rho_b[0])<deltam:
                        self.rhol=np.delete(self.rhol,1,0)
                        self.rhol=np.append(self.rhol,[self.rho_b[1]],0)
                    bestm=self.merit(qmax)
                    #print 'shifted'

                if bestm<bestbestm:
                    if np.abs(bestm-bestbestm)*100.0/bestbestm<0.1:
                        if repeat>1:
                            repeat=0
                            break
                        else:
                            repeat+=1
                    dm=bestbestm-bestm
                    bestbestm=bestm
                    bestrho=copy.copy(self.rhol)
                # For adjacent layer averaging to smooth out the randomness of the electron density
                tmprho=copy.copy(self.rhol)
                if avg>0 and np.mod(iterations,iprint)==0:
                    self.delta=self.delta*(1.0-5e-2)
                    for k in range(tnlayer,0,-1):
                        srho=0.0
                        nr=0
                        for kc in range(-int(avg/2),int(avg/2)+1):
                            if tnlayer+2>k+kc>-1 and kc!=0:
                                srho=srho+tmprho[k+kc]
                                nr=nr+1
                            elif k+kc<=-1:
                                srho=srho+self.rho_b[0]
                                nr=nr+1
                            elif tnlayer+2<=k+kc:
                                srho=srho+self.rho_b[1]
                                nr=nr+1
                        self.rhol[k]=srho/nr
                        bestm=self.merit(qmax)
                iterations=iterations+1
                if iprint!=0 and np.mod(iterations,iprint)==0:
                    print(iterations, bestm, bestbestm, self.delta,self.rhol[1])
                if iplot!=0 and np.mod(iterations,iplot)==0:
                    trho=copy.copy(self.rhol)
                    self.ax1.clear()
                    self.ax2.clear()
                    self.ax3.clear()
                    chisq=np.vstack((chisq,np.array([iterations,bestbestm])))
                    self.z=np.cumsum(self.dl)
                    #width=np.append(self.dl[1:],[0.0],0)
                    #self.z[-1]=self.z[-2]+self.dl[-2]
                    self.ax1.plot(self.x,self.y,'b.',label='Data')
                    self.ax1.set_xlabel('Q$_z$(\u212B$^{-1}$)')
                    self.ax1.set_ylabel('R/R$_F$')
                    self.ax1.set_title('Fresnel Normalized Reflectivity')
                    self.ax3.semilogy(chisq[:,0],chisq[:,1],'g*')
                    self.ax2.set_title('Electron Density Profile')
                    self.ax1.plot(self.q,self.func(self.q),'r-',label='Current fit')
                    self.ax2.step(self.z,(np.array(self.rhol)-self.rho_b[0])/(self.rho_b[1]-self.rho_b[0]),'r-',label='Current fit')
                    self.ax2.set_xlabel('Depth (\u212B)')
                    self.ax2.set_ylabel('$\\rho$/$\\rho_w$')
                    self.rhol=bestrho
                    self.ax1.plot(self.q,self.func(self.q),'g-',label='Best fit')
                    self.ax2.step(self.z,(np.array(self.rhol)-self.rho_b[0])/(self.rho_b[1]-self.rho_b[0]),'g-',label='Best fit')
                    self.ax1.set_yscale('log')
                    self.ax3.set_title('Fitting Info')
                    self.ax3.set_xlabel('Iterations')
                    self.ax3.set_ylabel('$\chi^2$')
                    self.ax1.legend(loc='best')
                    self.ax2.legend(loc='best')
                    pl.tight_layout()
                    pl.pause(0.001)
                    self.rhol=trho
                    if saveFig and dir is not None:
                        pl.savefig(dir+'Figure_%05d.png'%self.counter)
                        self.counter+=1
        except (KeyboardInterrupt, SystemExit):
            self.rhol=copy.copy(bestrho)


    def cal_t(self,qz,qc):
        tsc=np.abs(2*qz/(qz+ np.sqrt(qz**2-qc**2+0j)))**2
        return tsc
    
    def func(self,x):
        """
        Calculates the fitting function at the q values.
        x=list of q values
        """
        refq,r=parratt(x,self.lamda,self.dl,self.rhol,self.betal)
        rf,r=parratt(x,self.lamda,[0.0,0.0],self.rho_b,self.beta_b)
        if self.Rf!=0:
            refq=refq/rf
        if self.logy!=0:
            return np.log10(refq)
        else:
            return refq

    
    def phiqz(self,x):
        """
        Calculates structure factor for GIXOS calculation
        """
        if type(x)!=list:
            x=[x]
        z=np.cumsum(self.dl)
        drho=np.append([0],np.diff(self.rhol))
        phi=[]
        for q in x:
            phi.append(np.sum(drho*np.exp((0+1j)*q*z*self.dl[1])))
        return np.abs(np.array(phi))**2

    def capillary(self,qz,qpar,gam,T,qmax):
        """
        Calculates the diffuse scattering part of scattering due to capillary wave
        """
        KbT=1.38e-23*(273.0+T)
        eta=KbT*qz**2*1e23/2.0/np.pi/gam
        return (qpar*np.exp(-0.5772)/qmax)**eta*eta/qpar**2/qz**2

    def gixos(self,qz,qpar,gam,T,qmax):
        """
        Calculates the gixos as a funciton of qz at a particular qpar
        """
        gix=[]
        for qz1 in qz:
            gix.append(self.phiqz(qz1)[0]*self.capillary(qz1,qpar,gam,T,qmax)*self.cal_t(qz1,0.0217))
        return np.array(gix)

    def gixos_merit(self,qm,qpar,gam,T,qmax):
        """
        Computes the merit function for GIXOS fitting
        """
        x=self.x[np.argwhere(self.x<qm)[:,0]]
        gix=self.gixos(x,qpar,gam,T,qmax)
        return np.sum(((self.y[:len(x)]/self.y[0]-gix/gix[0])*self.y[0]/self.y_err[:len(x)])**2)/len(x)#*np.sum(np.sqrt(1.0+np.diff(self.rhol)**2))        

    def gixos_evolve(self,qpar,gam,T,qm=0.7,qmax=2.23,deltam=1,maxiter=10**6,iplot=0,iprint=0,avg=0):
        """
        Evolve the electron densities of the layers to fit the GIXOS data
        qpar=in-plane reciprocal wave-vector on which GIXOS is measured
        gam=Interfacial tension in mN/m
        T=Temperature in Degree Celcius
        qm = maximum q value of the data to be used for the fitting
        deltam=minimum change in merit function or chi square value for stopping the iterations
        iplot=0 or N for not or plotting the results after N iterations
        iprint=0 or N for not or printing the results after N iterations
        avg=Number of sublayers to be used for averaging 
         """
        bestm=self.gixos_merit(qm,qpar,gam,T,qmax)
        bestrho=copy.copy(self.rhol)
        bestd=copy.copy(self.dl)
        bestbeta=copy.copy(self.betal)
        dm=bestm
        iterations=0
        tnlayer=len(self.rhol)-2
        chisq=np.array([1,bestm])
        bestbestm=bestm
        bestrho=copy.copy(self.rhol)
        while iterations<maxiter and self.delta>deltam:
            try:  
                for j in range(1,tnlayer+1):
                    if tnlayer>1:
                        k=np.random.randint(1,tnlayer)
                    else:
                        k=1
                maxr,minr=self.minmax_rho(self.rhol, k, self.delta)
                if abs(maxr-minr)>1e-5:
                    np.random.seed()
                    temp_rho=minr+np.random.rand()*(maxr-minr)
                    if self.rho_range[1]>temp_rho>self.rho_range[0]:
                        old_rho=self.rhol[k]
                        self.rhol[k]=temp_rho
                        m1=self.gixos_merit(qm,qpar,gam,T,qmax)
                        if m1>=bestm:
                            self.rhol[k]=old_rho
                        else:
                            bestm=m1
                    while np.abs(self.rhol[1]-self.rho_b[0])<deltam:
                        self.rhol=np.delete(self.rhol,1,0)
                        self.rhol=np.append(self.rhol,[self.rho_b[1]],0)
                    bestm=self.gixos_merit(qm,qpar,gam,T,qmax)
                    #print 'shifted'

                if bestm<bestbestm:
                    if np.abs(bestm-bestbestm)*100.0/bestbestm<0.0001:
                        break
                    dm=bestbestm-bestm
                    bestbestm=bestm
                    bestrho=copy.copy(self.rhol)
                tmprho=copy.copy(self.rhol)
                if avg>0 and np.mod(iterations,5)==0:
                    self.delta=self.delta*(1.0-5e-2)
                    for k in range(tnlayer,0,-1):
                        srho=0.0
                        nr=0
                        for kc in range(-avg/2,avg/2+1):
                            if tnlayer+2>k+kc>-1 and kc!=0:
                                srho=srho+tmprho[k+kc]
                                nr=nr+1
                            elif k+kc<=-1:
                                srho=srho+self.rho_b[0]
                                nr=nr+1
                            elif tnlayer+2<=k+kc:
                                srho=srho+self.rho_b[1]
                                nr=nr+1
                        self.rhol[k]=srho/nr
                        bestm=self.merit(qmax)
                iterations=iterations+1
                if iprint!=0 and np.mod(iterations,iprint)==0:
                    print(iterations, bestm, bestbestm, self.delta,self.rhol[1])
                if iplot!=0 and np.mod(iterations,iplot)==0:
                    trho=copy.copy(self.rhol)
                    self.ax1.clear()
                    self.ax2.clear()
                    self.ax3.clear()
                    chisq=np.vstack((chisq,np.array([iterations,bestbestm])))
                    self.z=np.cumsum(self.dl)
                    #width=np.append(self.dl[1:],[0.0],0)
                    #self.z[-1]=self.z[-2]+self.dl[-2]
                    self.ax1.plot(self.x,self.y/np.max(self.y),'b.')
                    self.ax3.semilogy(chisq[:,0],chisq[:,1],'g*')
                    gix=self.gixos(self.q,qpar,gam,T,qmax)
                    self.ax1.plot(self.q,gix/np.max(gix),'r-')
                    self.ax2.plot(self.z,(np.array(self.rhol)-self.rho_b[0])/(self.rho_b[1]-self.rho_b[0]),'r-')
                    self.rhol=bestrho
                    gix=self.gixos(self.q,qpar,gam,T,qmax)
                    self.ax1.plot(self.q,gix/np.max(gix),'g-')
                    self.ax2.plot(self.z,(np.array(self.rhol)-self.rho_b[0])/(self.rho_b[1]-self.rho_b[0]),'g-')
                    self.fig.canvas.draw()
                    self.rhol=trho
            except (KeyboardInterrupt, SystemExit):
                self.rhol=copy.copy(bestrho)
                break
    
      
    def merit(self,qmax):
        """
        Calculates the merit function or Chi square upto the qvalue given by qmax
        """
        x=self.x[np.argwhere(self.x<qmax)[:,0]]
        return np.sum((self.y[:len(x)]-self.func(x))**2)#/self.y_err[:len(x)])**2)/len(x)#*np.sum(np.sqrt(1.0+np.diff(self.rhol)**2))
        

    
        
        
