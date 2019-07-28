
import numpy as np
import matplotlib.pyplot as plt

class CellGuy:
    def __init__(self, start_loc, j_init, beta, omega):
        self.position = start_loc
        self.j = j_init
        self.beta = beta
        self.omega = omega
        
        #these may be somewhere else - what is a best practice?
        self._length = 100
        self._rho_tot = 80
        self.end = self.position + self._length
        self.rhos = np.zeros(self._length,float)
        self.rhos[79:len(self.rhos)-1] = 0.5 * np.random.rand(20)
        
        self.contact = False #changes if in contact with another cell
        self.inh = np.zeros(self._length,float) #this can be updated upon contact
        self.diff_rho = 0.1
        self.diff_I = 0.1 #find bio relevant values
        #chem will be a variable in CellPlayground.py
        
    def react(self, rho_a, rho_b, chem):
        '''
        onrate=b*(k0+gamma*a**2/(K**2+a**2) - beta*inh + omega*chem)
        onrate[onrate<0]=0
        #this is the RDE, including an inhibitor that propels Rho away and a chemoattractant
        return (onrate - delta*a)
        '''
        delta=1     #
        gamma = 1   #  Initialize these above?
        K = 1       #  Would it make sense for different cells to have different
        k0=0.067    #     values for these?
        N = len(self.rhos)-1 
        onrate = rho_b * (k0 + gamma * rho_a**2 / (K**2 + rho_a**2) 
                - self.beta*self.inh[1:N] + self.omega*chem)
        onrate[onrate<0] = 0
        offrate = delta * rho_a
        return onrate - offrate
    
    
    def diffuse(self,chem):
        '''        Rhos[0]=(Rhos[2]+Rhos[1])/2  
        Rhos[N]=(Rhos[N-2]+Rhos[N-1])/2
        a=np.sum(Rhos)/L     
        b=(Ntot-a*L)/L

        Rhos2[1:N] = Rhos[1:N] + const*(Rhos[2:N+1]+Rhos[0:N-1]-2*Rhos[1:N]) + h*react_ab(Rhos[1:N],b)
        Rhos,Rhos2 = Rhos2,Rhos'''
        N = len(self.rhos)-1          #because Python indexing is weird
        rho_new = np.zeros(N+1,float) #
        asteps = self._length/N
        self.h=(asteps/(0.5*self.diff_rho)) *0.008
        const = self.h*self.diff_rho/asteps**2
        
        
        
        #Boundaries
        self.rhos[0] = (self.rhos[2]+self.rhos[1])/2
        self.rhos[N] = (self.rhos[N-1]+self.rhos[N-2])/2
        #Boundaries if in contact
        if self.contact != False:
            inh2 = np.zeros(N+1,float)
            #inward flux
            if self.contact == 'left': 
                self.inh[0] = self.j/self.diff_I + self.inh[1]   #(inh[2]+inh[1])/2
                self.inh[N]=(self.inh[N-1] + self.inh[N-2]) / 2
            elif self.contact == 'right':
                self.inh[0]=(self.inh[2]+self.inh[1]) / 2
                self.inh[N]=self.j/self.diff_I + self.inh[N-1]
            #diffuse inhibitor
            inh2[1:N] = self.inh[1:N] + const*(
                           self.inh[2:N+1]+self.inh[0:N-1]-2*self.inh[1:N])
            
            #inh2,self.inh = self.inh,inh2
            self.inh = inh2
            
            
        #conservation of Rho
        self.rho_a = np.sum(self.rhos)/self._length
        self.rho_b = (self._rho_tot - self.rho_a*self._length)/self._length
        
        #RDE
        rho_new[0] = self.rhos[0]; rho_new[N] = self.rhos[N]#changed from N-1
        rho_new[1:N] = self.rhos[1:N] + const*(self.rhos[2:N+1]+self.rhos[0:N-1]-2*self.rhos[1:N]) + self.h*self.react(self.rhos[1:N],self.rho_b,chem[self.position+1:self.end-1]) #might be too long?
        #rho_new += self.h*self.react(self.rhos[1:N],self.rho_b,chem[self.position:self.end])
        #rho_in,rho_new = rho_new,rho_in                                                h*react_CI(Rhos[1:N],b,inh[1:N],chem[1:N])
        self.rhos,rho_new = rho_new, self.rhos
        #return rho_new
    
            
    def move(self,rho_profile):
        self.velocity = rho_profile[len(rho_profile-1)] - rho_profile[0] 
        #This may not need to be an attribute - but it could be good for having
        #some sense of inertia, acceleration, force, etc.
        self.position += self.velocity #should be able to move left and right ?
        
        #Don't go past the environment bounds -- 
        #how do I do this for the right bound? 
        #Should I include the boundary somehow? Should the cell loop around?
        #(including this could mess with the gradient)
        
        if self.position < 0:
            self.position = 0
        self.end = self.position + self._length
        
        '''
        if self.end > len(environment) - 1:
            self.end = len(environment) - 1
            self.position = len(environment) - 1 - self._length
        '''
        
        
            
        
        
#Cell()

