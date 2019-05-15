# -*- coding: utf-8 -*-
"""
Created on Wed May 15 09:51:57 2019

@author: Mathilde

"""

import numpy as np  ### Never use * to import stuff, as it makes it difficult to retrace where functions come from
import scipy as sp
import math as m
from scipy import interpolate  ### Useful to interpolate stuff
from scipy import integrate
from matplotlib import pyplot as plt

from loading_and_moment_diagrams import load_diagrams, Loadcalculator
N = 100.
Mydistribution = load_diagrams(N)[2]
Mzdistribution = load_diagrams(N)[3]
Tdistribution = load_diagrams(N)[4]



def load_airfoil(filename):
    f = open(filename,'r')
    lines = f.readlines()
    data = []
    data_z = []
    data_y = []
    for line in lines:
        x = line.split()
        data.append(x)
    for i in range(len(data)):
        for j in range(len(data[i])):
            data[i][j] = float(data[i][j])
    for i in range(len(data)):
        data_z.append(data[i][0])
        data_y.append(data[i][1])

    #plt.plot(data_z,data_y)
    #plt.show()
    return data,data_z,data_y

#print(load_airfoil('NACA3414.txt')[1])

def wing_stress(Mydistribution, Mzdistribution, Tdistribution):
    I_xx = 1
    I_yy = 1
    I_zz = 0.9 
    I_xy = 1
    I_zy = 1
    I_xz = 1
    
    #loading the airfoil datapoints
    data = load_airfoil('NACA3414.txt') 
    data_z = data[1]
    data_y = data[2]
    b = 60.90
    alpha = 0.3
    N = 100 #discretization step in x direction
            # number of points evaluated in airfoil circumference
    HalfspanValues = np.linspace(0, b / 2 - 0.00001, N)
    # transform the moments to the airfoil reference system 
    My = []
    Mz = []
    Mx = []
    stress_x = []
    min_stress_x = []
    max_stress_x = []
    for i in range(len(HalfspanValues)):
        M_y = np.cos(alpha)*Mydistribution[i] + np.sin(alpha)*Mzdistribution[i]
        M_z = -np.cos(alpha)*Mzdistribution[i] + np.sin(alpha)*Mydistribution[i]
        My.append(M_y)
        Mz.append(M_z)
        for j in range(len(data_z)):
            stressx = ((M_z*I_yy - M_y*I_zy)*data_y[j] + (M_y*I_zz - M_z*I_zy)*data_z[j])/(I_zz*I_yy - I_zy**2)
            stress_x.append(stressx)
    stress_x = np.asarray(stress_x)
    stress_x = np.reshape(stress_x, (len(HalfspanValues),len(data_z)))

    for i in range(len(stress_x[:])):
        min_stress_x.append(min(stress_x[i]))
        max_stress_x.append(max(stress_x[i]))
        
    print(min_stress_x[0])
    print(max_stress_x[0])
    
    plt.subplot(2,1,2)
    plt.subplot(2,1,1)
    plt.plot(data_z,stress_x[0], 'r')
    plt.subplot(2,1,2)
    plt.plot(data_y,stress_x[0], 'g')
    plt.show()
        
        
print(wing_stress(Mydistribution, Mzdistribution, Tdistribution))    
    
    
    
    
