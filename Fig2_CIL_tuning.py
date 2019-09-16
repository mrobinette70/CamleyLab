# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 23:00:43 2019

@author: micha
"""

"""Notes:
    fix contact recognition: maybe just put it in another function
"""


import numpy as np
import matplotlib.pyplot as plt
from CellTest import CellTest
from CellFunctions import boundary_check, place_cells, build_timegraph, float_bound_check, check_contact
from scipy.stats import sem

'''This code models CIL in the presence of different chemoattractant
gradients (as in Lin et al Fig2). This should be wrapped in a function 
so it can be called in ParameterTuningMaster.
See Fig1c_w_v_tuning.py for advice on bizarre errors with the
try-except-finally block.
'''

env_length = 1500 #this isn't really spatially dependent since the gradient is uniform - cutting down on computational time with a smaller well
cell_env = np.zeros(env_length,float)
#env_chem = np.zeros(env_length,float)#
env_chem = np.linspace(2.5,5.8,env_length)  #nM
#env_chem = np.linspace(6.6,9.9,env_length)
#env_chem = np.linspace(10.8,14.1,env_length)

outcome1 = 0 #both cells move right (up gradient)
outcome2 = 0 #cells move apart
outcome3 = 0 #both cells move left
errors = 0   #just in case the cells don't separate - this should be 0 if the sim runs long enough


for i in range(1):
    print(i) #this takes a while to run - sanity check
    try:   
        c1 = CellTest(start_loc=300, polarity='right')
        c2 = CellTest(start_loc=1200, polarity='left')

    
        t = 0
        t_end = 10 #s
           
        cells = [c1,c2] #list of all initialized cells - make a function to make many?
    
        
        counter = 0
        increment = 5
        tg = build_timegraph(c1.h, t_end, env_length, every=increment)
        
        while t < t_end:
            for c in cells:
                print('t: ',t)
                c.diffuse(env_chem)
                c.move()
                
            t += c1.h
            float_bound_check(cells)
            check_contact(cells)

            
            if counter % increment == 0:
               ''' filling array to have cells over time '''
               allcells = place_cells(cells, env_length)
               tg[int(counter/increment),:] = allcells
                
            counter += 1
        
 
    except:
        continue
    finally: 
        if c1.left_tally == 0 and c2.right_tally > 0: #cells move right
            outcome1 += 1
        elif c1.left_tally > 0 and c2.right_tally > 0: #cells move apart
            outcome2 += 1
        elif c1.left_tally > 0 and c2.right_tally == 0: #cells move left
            outcome3 += 1
        elif c1.left_tally == 0 and c2.right_tally == 0: #some error
            errors += 1
        plt.imshow(tg)
print('Cells move up gradient (toward high concentration)',outcome1,'times.')
print('Cells move apart',outcome2,'times.')
print('Cells move down gradient',outcome3,'times.')
print('Cells do not separate',errors,'times.')