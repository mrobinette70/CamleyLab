import numpy as np
import matplotlib.pyplot as plt

class CellGuy:
    def __init__(self, start_loc, j_init, beta, omega):
        self.position = start_loc
        self.j = j_init
        self.beta = beta
        self.omega = omega
        
        self._length = 100
        self.end = self.position + self._length
        self.rhos = np.zeros(self._length,float)+.2 #how did I do this before?
        
    def diffuse(self):
        return
    
    def react(self,chem=False,contact=False):
        if contact == False:
            self.j=0
        elif contact == True:
            self.j=self.j_init
            
    def move(self,rho_profile):
        self.velocity = rho_profile[len(rho_profile-1)] - rho_profile[0]
        self.position += self.velocity
        
        #don't go past the environment bounds -- 
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

