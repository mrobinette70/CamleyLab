# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 15:06:54 2019

@author: micha
"""

import numpy as np
import matplotlib.pyplot as plt
from CellClass import CellGuy


env_length = 150
cell_env = np.zeros(env_length,float)
chem = np.zeros(env_length,float)
c1 = CellGuy(1,.1,.05,.05)
c2 = CellGuy(150,.1,.05,.05)

#plt.plot(c1.rhos,'.')

def vis_cells(cells,env_length):
    '''
    This function visualizes all cells in a given environment.
    
    cells: a tuple of all cells (from instance of CellClass)
    env_length: the length of the 1D environment cells swim in
    '''
    cell_env = np.zeros(env_length,float)
    #if len(cells) > 1:
     #   for cell in cells:  
      #      cell_env[cell.position:cell.end] = cell.rhos
    #else:
    cell_env[cells.position:cells.end] = cells.rhos
    plt.plot(cell_env,'.r')
    
#vis_cells((c1),150)
        

t = 0
t_end = 10

while t < t_end:
    #for cell in (c1):
     #  cell.diffuse(chem[cell.position:cell.end])
    c1.diffuse(chem)#[c1.position:c1.end])
    
    t += c1.h

vis_cells((c1),150)
        

'''
how can I include Rhos in the plot?

- have something like: cell_env[c1.position:c1.end] = c1.rhos for all cells,
then plot that outcome - maybe with different colors. This will maybe work
'''
'''
Problems:
    what's up with these plots???


'''

