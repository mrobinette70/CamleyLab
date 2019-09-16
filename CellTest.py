# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 19:11:57 2019

@author: micha

Things to add:
    Velocity calculations don't need to be rounded. Round when the position is needed
    to plot.
        This will also mean there needs to be a slight refactoring of boundary_check - maybe do it the 
        same way, but just do something like:
            if cell1 in cell2:
                cell1.position += (cell1.position - cell2.position)
    Velocity calculations need to be normalized by length in some way, so longer 
        cells aren't automatically faster. Maybe they should even be slower? Could divide 
        by length^2 as some function of drag with the environment.
    Add some noise to inhibitor influx - may look like something with h or h^1/2 variance, mean on 0 
    Change parameters to match timescale
    I probably need to make things in terms of min-1. For some reason the really small parameter
        values make this take forever to run.
    
                



"""

import numpy as np

class CellTest:
    
    #Shouldn't be changed - if there's a problem, one of these might be changed by mistake
    L=10                      # Cell length (um)
    Da=0.1                    # Rho diffusion constant (in membrane) (um2 s-1)
    Db=10                     # Rho diffusion constant (in cytosol) - not used in this implementation (approximation)
    Ntot=80                   # Total amount of Rho in cell (maybe this could be variable too?)
                              # Ntot was something smaller in the paper used - 22.6833 (change back before tuning parameters)

    
    N=100                     # Number of spatial steps
    asteps=L/N                # Spatial step
    h=(asteps/(0.5*Db)) *0.8  # Time step -- h <= a/(2D) 
    const=h*Da/asteps**2      # Constant for RDE (Rho diffusion)
    
    
    

    def __init__(self, start_loc, beta=0.07, 
                 j=.5, Di=10, omega=0.05, kchem = 4.15, kinh = 5, v_const = 1e-3, polarity='random'): #customizable inputs
        
        self.beta = beta
        self.j = 0
        self.j_contact = j #assuming the cell is not in contact, then if in contact this will be ~~ the influx of inh
        self.Di = Di
        self.omega = omega
        self.kchem = kchem
        self.kinh = kinh
        self.v_const = v_const
        
        self.position = start_loc
        self.end = self.position + self.N
        if polarity == 'random':
            self.Rhos = np.random.random(size=self.N+1)*.1 #to be changed to loading premade polarized cells (or as an option: left, right, random?)
        elif polarity == 'left':
            #noise can be added optionally to the initial profiles (pre-initialized)
            self.Rhos = np.load('left_polarized.npy')
            #self.Rhos += (np.random.random(len(self.Rhos))-0.5)*1.3 #adding a little randomness
            self.Rhos[self.Rhos < 0] = 0
        elif polarity == 'right':
            self.Rhos = np.load('right_polarized.npy')
            #self.Rhos += (np.random.random(len(self.Rhos))-0.5)*1.3 #adding a little randomness
            self.Rhos[self.Rhos < 0] = 0
        
        

        self.Rhos2 = np.zeros((self.N+1),float)
        self.inh=np.zeros((self.N+1),float)
        self.inh2=np.zeros((self.N+1),float)

        self.left_contact = False
        self.right_contact = False
        
        #Tallies (for parameter tuning)
        self.direction = None      #right or left
        self.new_direction = None  #right or left
        self.left_tally = 0        #if the cell changes to go left
        self.right_tally = 0       #if the cell changes to go right
        self.contact_tally = 0     #total number of cell contacts (may need to divide this by 2 at the end after summing all cells)
        
        cm = np.floor(self.N/2)
        dir_num = np.sum([(idx-cm)*val for idx,val in enumerate(self.Rhos)])
        if dir_num >= 0:
            self.direction = 'right'
        elif dir_num < 0:
            self.direction = 'left'
            
        self.positions = []  # Optionally used below in self.move()
        self.velocities = [] #
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
        k0=0.067 # s-1                 # Constants for RDE
        delta=1  #s-1                  # 
        gamma=1  #s-1                  # 
        K=1                            # 
        
        onrate = b * (k0 + gamma*a**2/(K**2 + a**2)
               - self.beta * (inh/(self.kinh + inh)) #this might be good enough - not sure how to model saturation without having data on it 
               + self.omega * (chem/(self.kchem+chem))) #adding saturation - based on empirical data for best gradient
        onrate[onrate < 0] = 0
        offrate = delta*a
        #this is the RDE, including an inhibitor that propels Rho away and a chemoattractant
        return (onrate - offrate)
    
    def diffuse(self,env_chem):
        ''' 
        Diffusing forward one iteration
        Parameters:
            
        env_chem: array-like
            Chemoattractant profile of the environment. The RDE uses the chemoattractant within the 
            length of the cell (chem) to determine chemotaxis
                      
        '''
        
        N = self.N #having this call self.N gets very messy and hard to read
        Rhos = self.Rhos
        Rhos2 = self.Rhos2
        inh = self.inh
        inh2 = self.inh2
        
        #Adding Gaussian noise
        inh += np.random.normal(loc=0,scale=np.sqrt(self.h),size=len(inh))
        inh[inh < 0 ] = 0
        
        #Grabbing an array of chemoattractant concentrations in the environment
        #(rounded since positions are no longer rounded)
        chem=env_chem[int(np.rint(self.position)):int(np.rint(self.end+1))]
        
        if self.left_contact == True or self.right_contact == True:
            self.j = self.j_contact
        else:
            self.j = 0 #checking if the cells detach - the influx will need to return to 0

       
        


        if self.right_contact == True & self.left_contact == True:
            inh[N]=self.j 
            inh2[0]=inh[0]; inh2[N]=inh[N]
            
            inh[0] = self.j/self.Di + inh[1] 
            inh[N] = self.j/self.Di + inh[N-1]
        elif self.left_contact == True:
            inh[0]=self.j 
            inh2[0]=inh[0]; inh2[N]=inh[N]            
            
            inh[0] = self.j/self.Di + inh[1]
            inh[N] = (inh[N-1] + inh[N-2])/2
        elif self.right_contact == True:
            inh[N]=self.j 
            inh2[0]=inh[0]; inh2[N]=inh[N]
            
            inh[0] = (inh[2] + inh[1])/2
            inh[N] = self.j/self.Di + inh[N-1]
        else:
            inh[0] = (inh[2] + inh[1])/2 
            inh[N] = (inh[N-1] + inh[N-2])/2 
            
        Rhos[0]=(Rhos[2]+Rhos[1])/2
        Rhos[N]=(Rhos[N-2]+Rhos[N-1])/2 
        
        #Conservation of mass
        a=np.sum(Rhos)/self.L          
        self.b=(self.Ntot-a*self.L)/self.L

        Rhos2[0]=Rhos[0]; Rhos2[N]=Rhos[N]
        #Solving RDE
        inh2[1:N] = inh[1:N] + self.const*(inh[2:N+1]+inh[0:N-1]-2*inh[1:N])
        Rhos2[1:N] = Rhos[1:N] + self.const*(Rhos[2:N+1]+Rhos[0:N-1]-2*Rhos[1:N]) + self.h*self.react_CI(Rhos[1:N],self.b,inh[1:N],chem[1:N])#[self.position:self.end-1])


        self.Rhos = Rhos2 #equivalent to a return     
        self.inh = inh2
        
    def move(self):
        '''
        This function calculates the velocity of a cell given its Rho profile.
        
        Problems:
            What if it reaches a wall? This throws an error now.
            Should a cell's velocity depend on a contacting cell? Maybe the velocity for 
            both cells will be the sum of both velocities until contact breaks.
            ^^ Should vel be calculated somewhere else then? Either as another method or 
            as a function in the environment that accesses rho?
        '''
        #Useful for monitoring a cell's position over time (with self.velocities below)
        #self.positions.append(self.position) 
        
        def calc_vel():

            cm = np.floor(self.N/2)
            self.vel = self.v_const*np.sum([(idx-cm)*val for idx,val in enumerate(self.Rhos)])
        
        calc_vel()
        
        #Useful for determining if a cell changes directions - comparing direction and new_direction
        if self.vel >= 0:
            self.new_direction = 'right'
        elif self.vel < 0:
            self.new_direction = 'left'
            
        if self.direction == 'left' and self.new_direction == 'right':
            self.right_tally += 1
        elif self.direction == 'right' and self.new_direction == 'left':
            self.left_tally += 1
            
        self.direction = self.new_direction
        self.new_direction = None
        
        #vel > 0: move right
        if self.vel > 0 and self.right_contact == True: #unless contacting another cell!
            self.vel = 0
        elif self.vel < 0 and self.left_contact == True: #same with moving left
            self.vel = 0
        else:
            self.position += self.vel
            self.end += self.vel
            
        #Use this to monitor a cell's velocity over time    
        #self.velocities.append(self.vel)