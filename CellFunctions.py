# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 19:53:24 2019

@author: micha
"""
import numpy as np

def boundary_check(cells, env_length):
    '''
    For each cell, checking if in contact with the other cells. If so, the cell 
    is in "contact"
    
    cells: a tuple of all cells (from instance of CellClass)
    env_length: the length of the 1D environment cells swim in 
    
    Start with the first cell and push other cells out
    
    Label by cell #, instead of all 1. Count all 2s in cell 1 and push cell 2 out by that length, etc.
    This doesn't account for cells being bordered by two cells - seems to be fine, but this doesn't explicitly account for 
    stacks of cells. Hasn't been a problem in simulations with only a few (<5) cells.
    '''
    cell_env = np.zeros(env_length)
    for idx,c in enumerate(cells):
        cell_env[c.position:c.end+1] = idx+1 #establishing all cell positions
    for c in cells:
        
        if cell_env[c.position-1] != 0: #there's a cell to the left
            cellnum = cell_env[c.end+1] #number of the other cell
            howmany = np.count_nonzero(cell_env[c.position:c.end+1] == cellnum)
        
            c.position += howmany #move me out of the other cell
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

    cells: a tuple of all cells (from instance of CellTest)
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