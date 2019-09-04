# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 15:06:54 2019

@author: micha
"""
"""
Questions:
    Should inh be wiped after contact ends? Or should the inhibitor be allowed to 
    keep diffusing?
    

"""


import numpy as np
import matplotlib.pyplot as plt
from CellTest import CellTest


env_length = 50000
cell_env = np.zeros(env_length,float)
#env_chem = np.zeros(env_length,float)#
env_chem = np.linspace(0,100,env_length)
c1 = CellTest(start_loc=7500,beta=.1,j=1,omega=.9)
c2 = CellTest(start_loc=500, beta=.1,j=1,omega=1)
c3 = CellTest(start_loc=12000,beta=.1,j=1,omega=0.01)
c4 = CellTest(start_loc=10000, beta=.1,j=1,omega=0.1)

c2.Rhos=np.zeros(c2.N+1)
c2.Rhos[81:c2.N+1] = 0.2
#plt.plot(c1.rhos,'.')

def boundary_check(cells, env_length):
    '''
    For each cell, checking if in contact with the other cells. If so, the cell 
    is in "contact"
    
    cells: a tuple of all cells (from instance of CellClass)
    env_
    
    Start with the first cell and push other cells out
    
    Label by cell #, instead of all 1. Count all 2s in cell 1 and push cell 2 out by that length
    This doesn't account for cells being bordered by two cells
    '''
    cell_env = np.zeros(env_length)
    for idx,c in enumerate(cells):
        cell_env[c.position:c.end+1] = idx+1 #establishing all cell positions
    for c in cells:
        
        if cell_env[c.position-1] != 0: #there's a cell to the left
            cellnum = cell_env[c.end+1] #number of the other cell
            howmany = np.count_nonzero(cell_env[c.position:c.end+1] == cellnum)
        
            c.position += howmany
            c.end += howmany
            c.left_contact = True
        else:
            c.left_contact = False
            
        if cell_env[c.end+1] != 0: #there's a cell to the right
            cellnum = cell_env[c.end+1] #number of the other cell
            howmany = np.count_nonzero(cell_env[c.position:c.end+1] == cellnum)
            c.position -= howmany
            c.end -= howmany
            c.right_contact = True
        else:
            c.right_contact = False
        
def place_cells(cells, env_length):
    '''
    This function takes all cell Rho profiles and collects them into a single array
    for plotting over time.

    cells: a tuple of all cells (from instance of CellClass)
    env_length: the length of the 1D environment cells swim in    
    '''
    cell_env = np.zeros(env_length,float)
    #if len(cells) > 1:
    for c in cells:  
        cell_env[c.position:c.end+1] = c.Rhos
    return cell_env

def build_timegraph(h, t, env_length, every=1):
    '''
    Builds the array used to show cells over time in a surface/heat map
    (time axis is axis 0)
    
    Parameters:
        
        h: float
            Timestep used in CellTest.diffuse(). Retrieve from an instance of 
            CellTest.
        t: float, int
            Time that will elapse over the period of the simulation
        env_length: float, int
            Length of the environment cells can swim in
        every: int
            How frequently the cell will be plotted. every=1 plots the cell at 
            every timestep.
    '''
    time_array = np.zeros([int(np.ceil(t/(h*every)))+1,env_length])
    return time_array
"""
def vis_cells(cells,env_length):
    '''
    This function visualizes all cells in a given environment.
    
    cells: a tuple of all cells (from instance of CellClass)
    env_length: the length of the 1D environment cells swim in
    '''

    plt.plot(cell_env,'.r')
"""    
#vis_cells((c1),150)
        

t = 0
t_end = 100
cells = [c1,c2,c3,c4]
#for cell in cells:

#plt.plot(c1.Rhos)

print('Start and end: ',c1.position,', ',c1.end)
print(' ')

counter = 0
tg = build_timegraph(c1.h, t_end, env_length, every=1)

while t < t_end:
    for c in cells:
    #cell.diffuse(chem[cell.position:cell.end])
        c.diffuse(env_chem)
        c.move()
        '''
        if c1.left_contact == True:
            plt.figure(1)
            plt.plot(c1.Rhos)
            plt.figure(2)
            plt.plot(c1.inh)
            print(c1.vel)
        if c2.right_contact == True:
            print(t)
        '''    
    #boundary_check([cell_list])
    t += c1.h
    boundary_check(cells,env_length)
    allcells = place_cells(cells, env_length)

    ''' filling array to have cells over time '''
    ###fill array at [counter,:]
    tg[counter,:] = allcells
    counter += 1
    

    #print('Velocity: ',c1.vel)
    #print('Start and end: ',c1.position,', ',c1.end)
    #print(' ')
        #print(cell.rho_a)# + cell.rho_b)
#plt.figure(2)
#vis_cells((c1),env_length)
#print(c1.rho_a+c1.rho_b)
#plt.figure(3)        
#plt.plot(c1.Rhos)
'''
how can I include Rhos in the plot?

- have something like: cell_env[c1.position:c1.end] = c1.rhos for all cells,
then plot that outcome - maybe with different colors. This will maybe work
'''
'''
Problems:
    what's up with these plots???


'''

