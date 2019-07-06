# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 15:06:54 2019

@author: micha
"""

import numpy as np
import matplotlib.pyplot as plt
from CellClass import CellGuy

x=CellGuy(1,2,3,4)._length

env_length = 1000
cell_env = np.zeros(env_length,float)
chem = np.zeros(len(cell_env),float)
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
    for _cell in cells:
        cell_env[_cell.position:_cell.end] = _cell.rhos
    plt.plot(cell_env,'.')
    
vis_cells((c1,c2),300)
        


'''
how can I include Rhos in the plot?

- have something like: cell_env[c1.position:c1.end] = c1.rhos for all cells,
then plot that outcome - maybe with different colors. This will maybe work
'''

