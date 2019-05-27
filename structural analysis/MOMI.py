# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:35:21 2019

@author: Mels
"""
#from loading_and_moment_diagrams import *
#from math import *
#from stress_distribution_wing import *
from airfoil_geometry import airfoil_geometry
## Moment of Inertia of a wing section





#% load the airfoil into the program:

airfoil = 'NACA3414.txt'

def s_airfoil(N,b):
    
    data_z_all_sec = airfoil_geometry(N,b)[0]
    data_y_upper_all_sec = airfoil_geometry(N,b)[1]
    data_y_lower_all_sec = airfoil_geometry(N,b)[2]
    ds_sec_all = []
    s_all_sec = []
    
    for i in range(len(data_z_all_sec)):
        for j in range(len(data_z_all_sec[0])-1):
            ds_sec = sqrt((data_z_all_sec[i][j+1]-data_z_all_sec[i][j])**2+(data_y_upper_all_sec[i][j+1]-data_y_upper_all_sec[i][j])**2) + sqrt((data_z_all_sec[i][j+1]-data_z_all_sec[i][j])**2+(data_y_lower_all_sec[i][j+1]-data_y_lower_all_sec[i][j])**2)
            ds_sec_all.append(ds_sec)
    
    ds_sec_all = np.asarray(ds_sec_all)
    ds_sec_all = np.reshape(ds_sec_all, (len(data_z_all_sec),len(data_z_all_sec[0])-1))
    
    for i in range(len(ds_sec_all)):
        s_all_sec.append(sum(ds_sec_all[i]))
        
        
    return(s_all_sec)














