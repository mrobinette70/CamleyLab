# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 15:06:54 2019

@author: micha
"""
"""
Questions:
    Should inh be wiped after contact ends? Or should the inhibitor be allowed to 
    keep diffusing?

Notes:
    In fourcells_good.png:
        Four cells, started from random Rho distributions. All four cells 
        (fortunately) oriented toward the center, so we get to see CIL for all 
        of them. No chemoattractant now, but testing tallying implementation in 
        CellTest class. 3 and 4 are working perfectly - c3 has 2 contacts, and 
        c4 has 1 (validated by inspecting the picture). Interestingly, c1 has 1155, 
        and c2 has 1156 - I'm assuming c1 would hit c2, direct slightly away, repolarize 
        and hit c2 again. These cells would probably be very close in vivo without touching 
        after the first separation, but the rounding done to calculate new positions
        pushed them back together. Left and right tallies also worked well, even for 
        c1 and c2. These numbers may be 1 higher than seen in the png, since the
        cells start at random distributions, and a cell may start looking like a left-facing
        cell due to randomness and polarize to the right.
        c1.left_tally = 1; c1.right_tally = 1
        c2.left_tally = 0; c2.right_tally = 1
        c3.left_tally = 1; c3.right_tally = 0
        c4.left_tally = 0; c4.right_tally = 1
        The +1 discrepancy is probably not an issue. This will be checked when starting
        from a polarized state - if it's still coming up then, the way I determine 
        direction or velocity may be wrong. Can also be checked by plotting Rhos 
        for all cells over time in one of these heatmaps, which is much easier to see.
        

"""


import numpy as np
import matplotlib.pyplot as plt
from CellTest import CellTest


env_length = 40000
cell_env = np.zeros(env_length,float)
env_chem = np.zeros(env_length,float)#
#env_chem = np.linspace(0,100,env_length)
c1 = CellTest(start_loc=7500,beta=.1,j=1,omega=.9)
c2 = CellTest(start_loc=10000, beta=.1,j=1,omega=0.1)
c3 = CellTest(start_loc=12000,beta=.1,j=1,omega=0.01)
c4 = CellTest(start_loc=25000, beta=.1,j=1,omega=1)



c1.Rhos=np.zeros(c2.N+1)
c1.Rhos[81:c2.N+1] = 0.2
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
            
            if c.left_contact == False: #the cell is being newly contacted
                c.contact_tally += 1
            c.left_contact = True
        else:
            c.left_contact = False
            
        if cell_env[c.end+1] != 0: #there's a cell to the right
            cellnum = cell_env[c.end+1] #number of the other cell
            howmany = np.count_nonzero(cell_env[c.position:c.end+1] == cellnum)
            
            c.position -= howmany
            c.end -= howmany
            
            if c.right_contact == False: #the cell is being newly contacted
                c.contact_tally += 1 
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
increment = 5
tg = build_timegraph(c1.h, t_end, env_length, every=increment)

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
    
    t += c1.h
    
        
    boundary_check(cells,env_length)
    
    if counter % increment == 0:
        ''' filling array to have cells over time '''
        allcells = place_cells(cells, env_length)
        tg[int(counter/increment),:] = allcells
        
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

