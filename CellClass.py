import numpy as np
import matplotlib.pyplot as plt

class CellGuy:
    def __init__(self, start_loc, j_init, beta, omega):
        self.position = start_loc
        self.j = j_init
        self.beta = beta
        self.omega = omega
        
        self._length = 100
        self._rho_tot = 80
        self.end = self.position + self._length
        self.rhos = np.zeros(self._length,float)+.2 #how did I do this before?
        
        self.contact = False #changes if in contact with another cell
        self.inh = np.zeros(self._length,float) #this can be updated upon contact
        self.diff_rho = 0.1
        self.diff_I = 0.1 #find bio relevant
        #chem will be a variable in CellPlayground.py
        
    def react(self, rho_a, rho_b, chem, inh):
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
    
   
        onrate = rho_b * (k0 + gamma * rho_a**2 / (K**2 + rho_a**2) 
                 - self.beta*inh + self.omega*chem)
        onrate[onrate<0] = 0
        offrate = delta * rho_a
        return onrate - offrate
    
    
    def diffuse(self,rho_in,inh,chem):
        '''        Rhos[0]=(Rhos[2]+Rhos[1])/2  
        Rhos[N]=(Rhos[N-2]+Rhos[N-1])/2
        a=np.sum(Rhos)/L     
        b=(Ntot-a*L)/L

        Rhos2[1:N] = Rhos[1:N] + const*(Rhos[2:N+1]+Rhos[0:N-1]-2*Rhos[1:N]) + h*react_ab(Rhos[1:N],b)
        Rhos,Rhos2 = Rhos2,Rhos'''
        len_ = len(rho_in)-1       #because Python indexing is weird
        rho_new = np.zeros(len_+1) #
        asteps = self._length/len_
        h=(asteps/(0.5*Db)) *0.8
        const = h*self.diff_rho/asteps**2
        
        
        
        #boundaries
        rho_in[0] = (rho_in[2]+rho_in[1])/2
        rho_in[len_] = (rho_in[len_-1]+rho_in[len_-2])/2
        if self.contact != False:
            inh2 = np.zeros(len_+1,float)
         #inward flux
            if self.contact == 'left': 
                inh[0] = self.j/self.diff_I + inh[1]   #(inh[2]+inh[1])/2
                inh[len_]=(inh[len_-1] + inh[len_-2]) / 2
            elif self.contact == 'right':
                inh[0]=(inh[2]+inh[1]) / 2
                inh[len_]=self.j/self.diff_I + inh[len_-1]
            inh2[1:len_] = inh[1:len_] + const*(
                           inh[2:len_+1]+inh[0:len_-1]-2*inh[1:len_])
            
            inh2,inh = inh,inh2
            
            
            
        #conservation of Rho
        rho_a = np.sum(rho_in)/self._length
        rho_b = (self._rho_tot - rho_a*self._length)/self._length
        
        #RDE
        rho_new[0] = rho_in[0]; rho_new[len_-1] = rho_in[len_-1]
        rho_new[1:len_] = rho_in[1:len_] + const*(
                          rho_in[2:len_+1]+rho_in[0:len_-1]-2*rho_in[1:len_]) + h*self.react(rho_a[1:len_],rho_b,inh[1:len_],chem[1:len_])
        rho_in,rho_new = rho_new,rho_in
        
        return rho_new
    
            
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

