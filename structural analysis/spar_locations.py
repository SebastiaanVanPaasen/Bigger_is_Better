# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 20:14:52 2019

@author: mathi
"""
import numpy as np
from loading_and_moment_diagrams import c
import loading_and_moment_diagrams as lm

N= 100
b = lm.b
HalfspanValues = np.linspace(0, b / 2 - 0.00001, N)

nr_spars = 4
first_spar = 0.2
last_spar = 0.6

def spar_loc(HalfspanValues, nr_spars, first_spar, last_spar):
    
    spar_loc_sec = []
    for i in range(len(HalfspanValues)):
        first_spar_location = first_spar*c(HalfspanValues[i])
        last_spar_location = last_spar*c(HalfspanValues[i])
        delta_spar = (last_spar_location-first_spar_location)/(nr_spars-1)
        spar_loc = np.arange(first_spar_location,last_spar_location+delta_spar,delta_spar)
        spar_loc_sec.append(spar_loc)
        
    return spar_loc_sec, delta_spar

spar_loc_sec, delta_spar = spar_loc(HalfspanValues, nr_spars, first_spar, last_spar)
