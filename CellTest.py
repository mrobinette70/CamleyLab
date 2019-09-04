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
                 j=0.05, Di=0.1, omega=0.7): #customizable inputs
        
        self.beta = beta
        self.j = 0
        self.j_contact = j #assuming the cell is not in contact, then if in contact this will be the influx of inh
        self.Di = Di
        self.omega = omega
        
        self.position = start_loc
        self.end = self.position + self.N
       
        self.Rhos = np.random.random(size=self.N+1)*.3
        #self.Rhos=np.zeros((self.N+1),float)
        #self.Rhos[81:self.N+1]=0.2
        #self.Rhos[0:20] = 0.2
        self.Rhos2 = np.zeros((self.N+1),float)
        self.inh=np.zeros((self.N+1),float)
        self.inh2=np.zeros((self.N+1),float)
        #self.contact = False #False, 'left', or 'right' or 'both'
        self.left_contact = False
        self.right_contact = False
        'maybe instead this should be left_contact or right_contact instead of contact'

    
    def react_CI(self, a, b, inh, chem):
        '''
        Reaction part of RDE - wave-pinning model from Mori et. al, balancing
        a positively-cooperative onrate and concentration-dependent linear offrate.
        Used in diffuse()
        
        Parameters:
            
        a: array-like
            Profile of active (membrane-bound) Rho proteins within the cell
        b: float
            Profile of inactive (cytosolic) Rho proteins within the cell - assumed
            to be constant in this approximation, since cytosolic proteins diffuse 
            several orders of magnitude faster than membrane-bound varieties.
        inh: array-like
            Profile of inhibitor molecules that cause contact-inhibition of locomotion
            (CIL). 
        chem: array-like
            Chemoattractant profile of the environment. The RDE uses the chemoattractant within the 
            length of the cell to determine chemotaxis
        '''    
        onrate=b*(self.k0+self.gamma*a**2/(self.K**2+a**2)
                - self.beta*inh + self.omega*chem)
        onrate[onrate<0]=0
        offrate = self.delta*a
        #this is the RDE, including an inhibitor that propels Rho away and a chemoattractant
        return (onrate - offrate)
    
    def diffuse(self,env_chem):
        ''' 
        Diffusing forward one iteration
        Parameters:
            
        env_chem: array-like
            Chemoattractant profile of the environment. The RDE uses the chemoattractant within the 
            length of the cell (chem) to determine chemotaxis
            
        To add:
            what happens at contact? It could just be
            self.jcontact = 0.05     in __init__(),
            self.j = 0 (no contact)
            if contact != False:
                j=self.jcontact
            This will make sure j can be customized, but it will start as 0 until
            contacting another cell (the details about 'left' and 'right' can be figured out later)
                
        '''
        
        N = self.N #having this call self.N gets very messy and hard to read
        Rhos = self.Rhos
        Rhos2 = self.Rhos2
        inh = self.inh
        inh2 = self.inh2
        
        
        chem=env_chem[self.position:self.end+1]
        
        if self.left_contact == True or self.right_contact == True:
            self.j = self.j_contact
        else:
            self.j = 0 #checking if the cells detach - the influx will need to return to 0

       
        Rhos2[0]=Rhos[0]; Rhos2[N]=Rhos[N]


        if self.right_contact == True & self.left_contact == True:
            inh[N]=self.j 
            inh2[0]=inh[0]; inh2[N]=inh[N]
            
            inh[0] = self.j/self.Di + inh[1] #check the sign here
            inh[N] = self.j/self.Di + inh[N-1]
        elif self.left_contact == True:
            inh[0]=self.j 
            inh2[0]=inh[0]; inh2[N]=inh[N]            
            
            inh[0] = self.j/self.Di + inh[1]
            inh[N] = inh[N-1]
        elif self.right_contact == True:
            inh[N]=self.j 
            inh2[0]=inh[0]; inh2[N]=inh[N]
            
            inh[0] = (inh[2] + inh[1])/2
            inh[N] = self.j/self.Di + inh[N-1]
        else:
            inh[0] = 0
            inh[N] = 0
            
        Rhos[0]=(Rhos[2]+Rhos[1])/2
        
        Rhos[N]=(Rhos[N-2]+Rhos[N-1])/2    
        a=np.sum(Rhos)/self.L     
        b=(self.Ntot-a*self.L)/self.L

        #diffusion and RD for inhibitor and Rho
        inh2[1:N] = inh[1:N] + self.const*(inh[2:N+1]+inh[0:N-1]-2*inh[1:N])
        #print('Chem len: ',len(chem[self.position:self.end-1]))
        Rhos2[1:N] = Rhos[1:N] + self.const*(Rhos[2:N+1]+Rhos[0:N-1]-2*Rhos[1:N]) + self.h*self.react_CI(Rhos[1:N],b,inh[1:N],chem[1:N])#[self.position:self.end-1])
 
        #Rhos,Rhos2 = Rhos2,Rhos
        #inh,inh2 = inh2,inh
        #Rhos[0] = Rhos[1]
        #Rhos[N] = Rhos[N-1]



        self.Rhos = Rhos2 #equivalent to a return     #, chem, inh
        self.inh = inh2   #profile of inhibitor
        
    def move(self):
        '''
        This function calculates the velocity of a cell given its Rho profile.
        Since the cell is at low Reynolds number, assume no momentum is kept over
        timesteps.
        
        Problems:
            What if it reaches a wall? This throws an error now.
            Should a cell's velocity depend on a contacting cell? Maybe the velocity for 
            both cells will be the sum of both velocities until contact breaks.
            ^^ Should vel be calculated somewhere else then? Either as another method or 
            as a function in the environment that accesses rho?
        '''
        self.vel = 100*(self.Rhos[self.N] - self.Rhos[0])
        self.roundvel = int(np.rint(self.vel))
        
        #vel > 0: move right
        if self.vel > 0 and self.right_contact == True: #unless contacting another cell!
            self.vel = 0
        elif self.vel < 0 and self.left_contact == True: #same with moving left
            self.vel = 0
        else:
            self.position += self.roundvel
            self.end += self.roundvel
            
    #def update(self):
     #   if                     
'''
t=0


while t < 100:

    c.diffuse(100)
    if int(t) % 20 == 0:
        plt.plot(c.Rhos)
    t += c.h
'''