# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 08:16:26 2019

@author: micha
"""
import numpy as np
import matplotlib.pyplot as plt
from CellTest import CellTest
from CellFunctions import boundary_check, float_bound_check, place_cells, build_timegraph
from scipy.stats import sem
from scipy.io import loadmat

def fig1c_check(omega_iter, v_const_iter, kchem_iter):

    env_length = 1500     #the well is ~150 um - this factor of 10 accounts for the 10 difference in cell length vs spatial steps
    #cell_env = np.zeros(env_length,float)
    #env_chem = np.zeros(env_length,float)
    env_chem = np.linspace(2.5,5.8,env_length)  #nM
    #env_chem = np.linspace(6.6,9.9,env_length)
    #env_chem = np.linspace(10.8,14.1,env_length)
    
    vels = np.zeros(100)        
    for i in range(100):
        c1 = CellTest(start_loc=1, omega=omega_iter, kchem=kchem_iter, v_const=v_const_iter, polarity='right')
        print(i)
        # ^ this is where the tuning happens
        #try: 
            
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
    
            vels[i] = c1.vel
            #except:
        #print('Position:', c1.position)
        
        
        #vels[i] = c1.vel
        #print('Error on: ',i)
        #continue
    """ Plotting over time - uncomment this guy    
    if counter % increment == 0:
        ''' filling array to have cells over time '''
        allcells = place_cells(cells, env_length)
        tg[int(counter/increment),:] = allcells
        
    counter += 1
    """

    print('Average vel: ', np.mean(vels[vels>0])) #there's a ~.5 chance to start going left, since the initial dist is random.
    print('std: ',np.std(vels[vels>0]))
    print('S.E.M.: ',sem(vels[vels>0])) #if it goes left, it will throw an error and have a negative velocity, but a cell moving 
                                        #right should never have a negative velocity - this grabs only the velocities needed
    print('Number viable: ',len([x for x in vels if x > 0]))