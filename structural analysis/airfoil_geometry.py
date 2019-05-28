# -*- coding: utf-8 -*-
"""
Created on Thu May 16 09:55:34 2019

@author: Mathilde
"""

import numpy as np  ### Never use * to import stuff, as it makes it difficult to retrace where functions come from
import scipy as sp
import math as m
from scipy import interpolate ### Useful to interpolate stuff
from scipy import integrate
from matplotlib import pyplot as plt
from stress_distribution_wing import load_airfoil
from loading_and_moment_diagrams import c

N = 100 
b = 60.#39.56 #41.76


def airfoil_geometry(N,b):

    HalfspanValues = np.linspace(0, b / 2 - 0.00001, N)
#    print(HalfspanValues[0])
    
    data_z_all_sec = []
    data_y_lower_all_sec = []
    data_y_upper_all_sec = []
    
    for i in range(len(HalfspanValues)):
        
        data_z, data_y = load_airfoil('NACA3414.txt')[1], load_airfoil('NACA3414.txt')[2] 
        data_z_order =  np.array(data_z[0:int((len(data_y)/2))+1])*c(HalfspanValues[i])
        data_z_all_sec.append(data_z_order)
        data_y_lower = np.array(data_y[(int((len(data_y)/2))):][::-1])*c(HalfspanValues[i])
        data_y_lower_all_sec.append(data_y_lower)
        data_y_upper = np.array(data_y[0:int((len(data_y)/2))+1])*c(HalfspanValues[i])
        data_y_upper_all_sec.append(data_y_upper)
    
    data_z_all_sec = np.asarray(data_z_all_sec)
    data_z_all_sec = np.reshape(data_z_all_sec, (len(HalfspanValues),len(data_z_order)))
    data_y_upper_all_sec = np.asarray(data_y_upper_all_sec)
    data_y_upper_all_sec = np.reshape(data_y_upper_all_sec, (len(HalfspanValues),len(data_z_order)))
    data_y_lower_all_sec = np.asarray(data_y_lower_all_sec)
    data_y_lower_all_sec = np.reshape(data_y_lower_all_sec,(len(HalfspanValues),len(data_z_order)))
    
#    print(data_y_upper_all_sec[0])
#    print(data_y_lower_all_sec[10])
    return data_z_all_sec, data_y_upper_all_sec, data_y_lower_all_sec


#data = load_airfoil('NACA3414.txt')
data_z_all_sec, data_y_upper_all_sec, data_y_lower_all_sec = airfoil_geometry(N,b)
        