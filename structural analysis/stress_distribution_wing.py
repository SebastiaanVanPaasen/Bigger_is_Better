# -*- coding: utf-8 -*-
"""
Created on Wed May 15 09:51:57 2019

@author: Mathilde

"""

import numpy as np  ### Never use * to import stuff, as it makes it difficult to retrace where functions come from

#from loading_and_moment_diagrams import load_diagrams

#Mydistribution = load_diagrams()[2]
#Mzdistribution = load_diagrams()[3]
#Tdistribution = load_diagrams()[4]



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


def wing_stress(Mydistribution, Mzdistribution, Tdistribution):
    I_xx =  1e-9
    I_yy =  1e-9
    I_zz = 1e-9
    I_xy = 1e-9
    I_zy = 1e-9
    I_xz = 1e-9
    
    alpha = 0.3
    N = 100 #discretization step in x direction
    nodes = 100 # number of points evaluated in airfoil circumference
    HalfspanValues = np.linspace(0, b / 2 - 0.00001, N)
    # transform the moments to the airfoil reference system 
    My = []
    Mz = []
    Mx = []
    for i in range(len(HalfspanValues)):
        My = np.cos(alpha)*Mydistribution[i] + np.sin(alpha)*Mzdistribution[i]
        Mz = -np.cos(alpha)*Mzdistribution[i] + np.sin(alpha)*Mydistribution[i]
        for j in range(len(nodes)):
            stress_x = ((Mz*Iyy - My*Izy)*y[j] + (My*Izz - Mz*Izy)*z[j])/(Izz*Iyy - Izy**2)
            
            
        
    
    
    
    
    
