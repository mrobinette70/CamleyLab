# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 19:11:57 2019

@author: micha
"""

import numpy as np

class CellTest:
    #initialize Rhos as a random dist (or something)
    
    #Shouldn't be changed - if there's a problem, one of these might be changed by mistake
    L=10                      # Cell length
    Da=0.1                    # Rho diffusion constant (in membrane)
    Db=10                     # Rho diffusion constant (in cytosol) - not used in this implementation
    Ntot=80                   # Total amount of Rho in cell (maybe this could be variable too?)
                              # Ntot was something smaller in the paper used - ~28 (change back before tuning parameters)
    
    k0=0.067                  # Constants for RDE
    delta=1                   # 
    gamma=1                   # 
    K=1                       # 
    
    N=100                     # Number of spatial steps
    asteps=L/N                # Spatial step
    h=(asteps/(0.5*Db)) *0.8  # Time step -- h <= a/(2D) 
    const=h*Da/asteps**2      # Constant for RDE (Rho diffusion)
    #const_inh=h*Di/asteps**2  # Constant for RDE (inhib diffusion) - should based on smaller number (Da here, may change)
    
    
     # Creates a random distribution of Rho in cell
    
    #def rho_CI(Rhos_in,beta=0.05,j=0.05,Di=0.1,omega=0.7,chem_max=0.5): #all of the parameters I'd like to test
    def __init__(self, start_loc, beta=0.05, 
                 j=0.05, Di=0.1, omega=0.7, 
                 chem_max=0.5): #customizable inputs
        
        self.beta = beta
        self.j = j
        self.Di = Di
        self.omega = omega
        self.chem_max = 0.5 #this should be in the cell's environment 
        self.position = start_loc
        
        self.end = start_loc + self.L
        
        #self.Rhos = np.random.random(size=self.N+1)
        self.Rhos=np.zeros((self.N+1),float)
        self.Rhos[81:self.N+1]=0.2
        self.Rhos2 = np.zeros((self.N+1),float)
        self.inh=np.zeros((self.N+1),float)
        self.inh2=np.zeros((self.N+1),float)
        self.contact = False #False, 'left', or 'right'

    
    def react_CI(self, a, b, inh, chem):
        onrate=b*(self.k0+self.gamma*a**2/(self.K**2+a**2)
                - self.beta*inh + self.omega*chem)
        onrate[onrate<0]=0
        offrate = self.delta*a
        #this is the RDE, including an inhibitor that propels Rho away and a chemoattractant
        return (onrate - offrate)
    
    def diffuse(self,tEnd):
        # Diffusing forward one iteration
        
        N = self.N #having this call self.N gets very messy and hard to read
        Rhos = self.Rhos
        Rhos2 = self.Rhos2
        inh = self.inh
        inh2 = self.inh2
        """introduction of inhibitor and chemoattractant"""
        #Rhos2=np.zeros((N+1),float)
        #inh=np.zeros((N+1),float)
        #inh2=np.zeros((N+1),float)
        chem=np.zeros((N+1),float)
        
        #chem=np.linspace(chem_max,0,101) 
        ####chem=np.linspace(0,self.chem_max,101)
       
        Rhos2[0]=Rhos[0]; Rhos2[N]=Rhos[N]
        inh[N]=self.j #some arbitrary concentration to start with
        inh2[0]=0; inh2[N]=inh[N]
        #t=0
        #tEnd2=200 #change this to 200
    
#        while t<tEnd: # this needs to be pulled out - the while loop will be called in the environment
#            Rhos[0] = (Rhos[2] + Rhos[1])/2    
#            Rhos[N] = (Rhos[N-2] + Rhos[N-1])/2
#            #Rhos[Rhos<0]=0 #avoids negative values -- all active Rho (a) is gone
#            inh[0] = (inh[2] + inh[1])/2
#            inh[N] = self.j/self.Di + inh[N-1]
#            
#            a=np.sum(self.Rhos)/self.L     
#            b=(self.Ntot-a*self.L)/self.L
#            
#            #diffusion and RD for inhibitor and Rho
#            inh2[1:N] = inh[1:N] + self.const*(inh[2:N+1]+inh[0:N-1]-2*inh[1:N])
#            Rhos2[1:N] = Rhos[1:N] + self.const*(Rhos[2:N+1]+Rhos[0:N-1]-2*Rhos[1:N]) + self.h*self.react_CI(Rhos[1:N],b,inh[1:N],chem[1:N])
#            #Rhos2[Rhos2<0]=0
#            Rhos,Rhos2 = Rhos2,Rhos
#            inh,inh2 = inh2,inh
#            t += self.h
        
        #Rhos[0] = (Rhos[2] + Rhos[1])/2    
        #Rhos[N] = (Rhos[N-2] + Rhos[N-1])/2
      
        #Rhos[0] = Rhos[1]
        #Rhos[N] = Rhos[N-1]
        
        inh[0] = (inh[2] + inh[1])/2
        inh[N] = self.j/self.Di + inh[N-1]
        Rhos[0]=(Rhos[2]+Rhos[1])/2  
        Rhos[N]=(Rhos[N-2]+Rhos[N-1])/2    
        a=np.sum(Rhos)/self.L     
        b=(self.Ntot-a*self.L)/self.L
            
        #diffusion and RD for inhibitor and Rho
        inh2[1:N] = inh[1:N] + self.const*(inh[2:N+1]+inh[0:N-1]-2*inh[1:N])
        Rhos2[1:N] = Rhos[1:N] + self.const*(Rhos[2:N+1]+Rhos[0:N-1]-2*Rhos[1:N]) + self.h*self.react_CI(Rhos[1:N],b,inh[1:N],chem[1:N])
 
        #Rhos,Rhos2 = Rhos2,Rhos
        #inh,inh2 = inh2,inh        
        #Rhos[0] = Rhos[1]
        #Rhos[N] = Rhos[N-1]
        

        
        self.Rhos = Rhos2 #equivalent to a return     #, chem, inh
        self.inh = inh2   #profile of inhibitor