# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 15:06:54 2019

@author: micha
"""
#keeping these here for my sanity
import numpy as np
import matplotlib.pyplot as plt
from CellTest import CellTest
from CellFunctions import boundary_check, float_bound_check, place_cells, build_timegraph
from scipy.stats import sem
from scipy.io import loadmat


'''
All of these functions simulate a single cell moving in a 150 um well. Used to model
cell swimming alone, for tuning omega, v_const, and kchem

The functions contain try-except blocks because the cells would sometimes swim
past either end and cause an error. This shouldn't matter for the sake of the simulation, so 
exceptions are passed to the "except", allowing the velocity to be collected. If there
are other errors, it will keep saying the cell has some weird velocity and finish very
quickly. In that case, remove the try-except and try again. This will allow 
errors to be seen and fixed.
'''


def fig1c_low(omega_iter, v_const_iter, kchem_iter):

    env_length = 1500     #the well is ~150 um - this factor of 10 accounts for the 10 difference in cell length vs spatial steps (L vs N in CellTest)
    env_chem = np.linspace(2.5,5.8,env_length)  #nM

    c1 = CellTest(start_loc=1, omega=omega_iter, kchem=kchem_iter, v_const=v_const_iter, polarity='right')
    # ^ this is where the tuning happens
    try: 
        
        t = 0
        t_end = 120 
    
        cells = [c1]
    
        """ Plotting over time - uncomment this guy 
        counter = 0
        increment = 5
        tg = build_timegraph(c1.h, t_end, env_length, every=increment)
        """
        while t < t_end:
            for c in cells:
            
                c.diffuse(env_chem)
                c.move() 
            
            t += c1.h
            
            float_bound_check(cells)
            #vels[i] = c1.vel #used in the print below - need to remake the vels array if you want to use this
        return c1.vel
            
    except:
 
        return c1.vel
        

    """ Plotting over time - uncomment this guy    
    if counter % increment == 0:
        ''' filling array to have cells over time '''
        allcells = place_cells(cells, env_length)
        tg[int(counter/increment),:] = allcells
        
    counter += 1
    """
'''These can be moved somewhere else later 
print('Average vel: ', np.mean(vels[vels>0])) #there's a ~.5 chance to start going left, since the initial dist is random.
print('S.E.M.: ',sem(vels[vels>0])) #if it goes left, it will throw an error and have a negative velocity, but a cell moving 
                                    #right shouldn't have a negative velocity - this grabs only the velocities needed
print('Number viable: ',len([x for x in vels if x > 0]))
'''

def fig1c_med(omega_iter, v_const_iter, kchem_iter):

    env_length = 1500     #the well is ~150 um - this factor of 10 accounts for the 10 difference in cell length vs spatial steps
    env_chem = np.linspace(6.6,9.9,env_length)

    c1 = CellTest(start_loc=1, omega=omega_iter, kchem=kchem_iter, v_const=v_const_iter, polarity='right')
    # ^ this is where the tuning happens
    try: 
        
        t = 0
        t_end = 120 
    
        cells = [c1]
    
        """ Plotting over time - uncomment this guy
        counter = 0
        increment = 5
        tg = build_timegraph(c1.h, t_end, env_length, every=increment)
        """
        while t < t_end:
            for c in cells:
            
                c.diffuse(env_chem)
                c.move() 
            
            t += c1.h
            
            float_bound_check(cells)
        return c1.vel
            #vels[i] = c1.vel
    except:
        return c1.vel
    
    
def fig1c_high(omega_iter, v_const_iter, kchem_iter):

    env_length = 1500     #the well is ~150 um - this factor of 10 accounts for the 10 difference in cell length vs spatial steps
    env_chem = np.linspace(10.8,14.1,env_length)

    c1 = CellTest(start_loc=1, omega=omega_iter, kchem=kchem_iter, v_const=v_const_iter, polarity='right')
    # ^ this is where the tuning happens
    try: 
        
        t = 0
        t_end = 120 
    
        cells = [c1]
    
        """ Plotting over time - uncomment this guy
        counter = 0
        increment = 5
        tg = build_timegraph(c1.h, t_end, env_length, every=increment)
        """
        while t < t_end:
            for c in cells:
            
                c.diffuse(env_chem)
                c.move() 
            
            t += c1.h
            
            float_bound_check(cells)
        return c1.vel
            #vels[i] = c1.vel
    except:
        return c1.vel    
    




