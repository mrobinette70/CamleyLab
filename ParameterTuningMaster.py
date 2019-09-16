# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 20:29:09 2019

@author: micha
"""

import numpy as np
import matplotlib.pyplot as plt
from CellTest import CellTest
from CellFunctions import boundary_check, float_bound_check, place_cells, build_timegraph
from scipy.stats import sem
from scipy.io import loadmat
from Fig1c_w_v_tuning import fig1c_low, fig1c_med, fig1c_high
from scipy.optimize import least_squares, minimize
from tqdm import tqdm
import datetime
'''
Ideal values for experiment 1:
    v(2.5-5.8 nM) = .016 um/s
    v(6.6-9.9 nM) = .015 um/s
    v(10.8-14.1 nM) = .0125 um/s
Difference from these values must be minimized
'''
print(datetime.datetime.now())
def chem_only(x):
    om = x[0] #omega
    vc = x[1] #v_const (the coefficient multiplied in velocity calculation)
    kc = x[2] #k_chem  (the term in the chemoattractant saturation of react_ab)
    expected = np.array([.016,.015,.0125])
    low = fig1c_low(om,vc,kc)
    med = fig1c_med(om,vc,kc)
    high = fig1c_high(om,vc,kc)
    
    res = (expected - np.array([low,med,high]))
    #print(np.array([low,med,high]))
    return res.dot(res)

x0chem = [2.47843255e-05, 4.20822806e-06, 3.54670223e-05]
    #x: array([2.47373751e-05, 4.23110681e-06, 3.53876513e-05]) #the output from the second experiment
out = minimize(chem_only,x0chem)






print(datetime.datetime.now()) #this can take a while - sanity check