#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 15:23:49 2019

@author: Max


"""

import numpy as np
import matplotlib.pyplot as plt
N_stringers = 100 #number of stringers for fuselage
r = 2.5 #fuselage diameter in meter
t_skin = 1. #skinn tickness in mm


stringer = []

for i in range(N_stringers):
    #make stringer list with a list of 3 parameters per stringer
    # stringer[i] has a list: [stringer_number, x_loc, y_loc]
    angle = np.deg2rad(360/N_stringers)
    x_loc = r*np.cos(i*angle)
    y_loc = r*np.sin(i*angle)
    stringer.append([i,x_loc,y_loc])
    
for i in range(len(stringer)):
    plt.plot(stringer[i][1],stringer[i][2],'x', label = str(stringer[i][0]))
plt.legend()
plt.show()    

