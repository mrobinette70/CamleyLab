# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 15:06:54 2019

@author: micha
"""
"""
Questions:
    Should inh be wiped after contact ends? Or should the inhibitor be allowed to 
    keep diffusing?
    Is there a better way to tune parameters than just guessing? With 3 separate 
    gradients, filling two parameters should be doable analytically, but the equations
    aren't a simple linear system.

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
from CellFunctions import boundary_check, place_cells, build_timegraph
from scipy.stats import sem

env_length = 1500     #the well is ~150 um - this factor of 10 accounts for the 10 difference in cell length vs spatial steps
cell_env = np.zeros(env_length,float)
#env_chem = np.zeros(env_length,float)#
#env_chem = np.linspace(2.5,5.8,env_length)  #nM
#env_chem = np.linspace(6.6,9.9,env_length)
env_chem = np.linspace(10.8,14.1,env_length)

#c1 = CellTest(start_loc=1,beta=.1,j=1)

#c1.Rhos=np.zeros(c1.N+1)
#c1.Rhos[0:20] = 0.2
#plt.plot(c1.rhos,'.')


vels = np.zeros(500)        
for i in range(500):
    c1 = CellTest(start_loc=1,beta=.1,j=1,polarity='right')
    try: 
        
        t = 0
        t_end = 120 * c1.h

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
            
            boundary_check(cells,env_length)
        vels[i] = c1.vel
    except:
        vels[i] = c1.vel
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
