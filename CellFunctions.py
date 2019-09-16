# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 19:53:24 2019

@author: micha
"""
import numpy as np

def float_bound_check(cells):
    '''
    This function checks if cells overlap, but it allows cell positions to be floats, making
    velocity calculations less prone to error from rounding.
    The function may be poorly optimized for large amounts of cells with O(n**2) but
    it should be fine for small groups used here.
    
    Interesting problem I found: The cells would overlap, and the way this works the 
    first cell would be pushed out over and over again, so the left cell would get pushed off the 
    boundary without actually having a negative velocity!  OR (and?) the first part needs to be there
    for this to work really well. Fix immediately (FIXED)
    
    This depends on which cell is first -- very bad (FIXED)
    '''
    for current_cell in cells:
        for other_cell in cells:

            if current_cell == other_cell: 
                pass
            elif (current_cell.end > other_cell.position) and (current_cell.end < other_cell.end):
                #the cell should go left
                #cells are both pushed outward equally
                overlap = (current_cell.end - other_cell.position)
                current_cell.position -= overlap/2
                current_cell.end = current_cell.position + current_cell.N
                other_cell.position += overlap/2 
                other_cell.end = other_cell.position + other_cell.N
            else:
                #cells aren't contacting
                pass
            
def check_contact(cells):
    '''
    This function replaces the part of float_bound_check that does contact 
    checking due to some bugs I couldn't pin down. Checks if there are any cells with 
    other cells at the same boundary positions and changes cell contact attributes.
    
    Care that eps is sufficiently small! If velocities are too low, having a small-ish eps might
    not be good enough - it should be something like 1e-2*velocity at least
    
    Can I do velocity calculations here, if the cells should have some net nonzero
    velocity after contact?
    '''
    for current_cell in cells:
        for other_cell in cells:
            eps = 1e-8 #probably small enough
            if current_cell == other_cell:
                #same cell
                pass
            elif np.abs(current_cell.position - other_cell.end) < eps:
                #other cell on the left of the current cell
                current_cell.left_contact = True
                other_cell.right_contact = True
            elif np.abs(current_cell.position - other_cell.end) > eps:
                current_cell.left_contact = False
                other_cell.right_contact = False
            
            if current_cell == other_cell:
                #same cell
                pass
            elif np.abs(current_cell.end - other_cell.position) < eps:
                #other cell on the right of the current cell
                current_cell.right_contact = True
                other_cell.left_contact = True    
            elif np.abs(current_cell.end - other_cell.position) > eps:
                #other cell on the right of the current cell
                current_cell.right_contact = False
                other_cell.left_contact = False                
                

#For graphing
        
def place_cells(cells, env_length):
    '''
    This function takes all cell Rho profiles and collects them into a single array
    for plotting over time.

    cells: a tuple of all cells (from instance of CellTest)
    env_length: the length of the 1D environment cells swim in    
    '''
    cell_env = np.zeros(env_length,float)
 
    for c in cells:  
        cell_env[int(np.rint(c.position)):int(np.rint(c.end+1))] = c.Rhos
    return cell_env


def build_timegraph(h, t, env_length, every=1):
    '''
    Builds the array used to show cells over time in a surface/heat map
    (time axis is axis 0)
    
    Parameters:
        
        h: float
            Timestep used in CellTest.diffuse(). Retrieve from an instance of 
            CellTest (c1.h, for example).
        t: float, int
            Time that will elapse over the period of the simulation
        env_length: float, int
            Length of the environment cells can swim in
        every: int
            How frequently the cell will be plotted. every=1 plots the cell at 
            every timestep (may be initialized as 'increment' elsewhere this is used).
    '''
    time_array = np.zeros([int(np.ceil(t/(h*every)))+1,env_length])
    return time_array



def boundary_check(cells, env_length): #no longer used - see float_bound_check 
    '''
    For each cell, checking if in contact with the other cells. If so, the cell 
    is in "contact"
    
    cells: a tuple of all cells (from instance of CellClass)
    env_length: the length of the 1D environment cells swim in 
    
    Start with the first cell and push other cells out
    
    Label by cell #, instead of all 1. Count all 2s in cell 1 and push cell 2 out by that length, etc.
    This doesn't account for cells being bordered by two cells - seems to be fine, but this doesn't explicitly account for 
    stacks of cells. Hasn't been a problem in simulations with only a few (<5) cells.
    
    ***No longer used, since position is no longer rounded. Kept for reference***
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
                
    
    
    
    
    
    
    
    
    
    