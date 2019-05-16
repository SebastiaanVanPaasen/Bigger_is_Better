# -*- coding: utf-8 -*-
"""
Created on Wed May 15 16:01:47 2019

@author: floyd
"""
import numpy as np
from stress_distribution_wing import load_airfoil


data = load_airfoil("naca3414.txt")
dataz = np.asarray(data[1])
datay = np.asarray(data[2])
Cr = 12.227

def y_upper(z):    
    y_upper = .5-z**2 #put in right y function
    return(y_upper)
def y_lower(z):
    y_lower = -.5+z**2 #put in right y function
    return(y_lower)

nr_spars = 4
first_spar_location = 0.2
last_spar_location = 0.75
delta_spar = (last_spar_location-first_spar_location)/(nr_spars-1)
spar_loc = np.arange(first_spar_location,last_spar_location+delta_spar,delta_spar)

upper_y = np.zeros(len(spar_loc))
lower_y = np.zeros(len(spar_loc))

for i in range(len(spar_loc)):
    upper_y[i] = y_upper(spar_loc[i])
    lower_y[i] = y_lower(spar_loc[i])
delta_y = abs(upper_y) + abs(lower_y)

spar_t = 0.002 #mm
I_zz_spar = 1/12*spar_t*delta_y**2 #+ steiner term, need centroid
I_yy_spar = 1/12*spar_t**2*delta_y #+steiner term
    
