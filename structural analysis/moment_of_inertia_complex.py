# -*- coding: utf-8 -*-
"""
Created on Wed May 15 16:01:50 2019

@author: Mathilde
"""
import numpy as np  ### Never use * to import stuff, as it makes it difficult to retrace where functions come from
import scipy as sp
import math as m
from stress_distribution_wing import load_airfoil
from loading_and_moment_diagrams import c

N = 100 
b = 60.90
HalfspanValues = np.linspace(0, b / 2 - 0.00001, N)
print(HalfspanValues[0])

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
#print(data_z_order)
#print(data_y_upper)
#print(data_y_lower)
data_z_all_sec = np.asarray(data_z_all_sec)
data_z_all_sec = np.reshape(data_z_all_sec, (len(HalfspanValues),len(data_z_order)))
data_y_upper_all_sec = np.asarray(data_y_upper_all_sec)
data_y_upper_all_sec = np.reshape(data_y_upper_all_sec, (len(HalfspanValues),len(data_z_order)))
data_y_lower_all_sec = np.asarray(data_y_lower_all_sec)
data_y_lower_all_sec = np.reshape(data_y_lower_all_sec,(len(HalfspanValues),len(data_z_order)))

#print(data_y_lower_all_sec)

t = 0.005

#print(c(HalfspanValues[0]))

def airfoil_geometry():
    
    for i in range(len(HalfspanValues)):
        Polyfit_airfoil_upper = sp.interpolate.interp1d(data_z_all_sec[i], data_y_upper_all_sec[i], kind="cubic", fill_value="extrapolate")
        Polyfit_airfoil_lower = sp.interpolate.interp1d(data_z_all_sec[i], data_y_lower_all_sec[i], kind="cubic", fill_value="extrapolate")
        
        N = 100.
        dz = c(HalfspanValues[i])/N
        #print(c(HalfspanValues[i]))
        #print(dz)
        z=0.
        start_wingbox = 0.2*c(HalfspanValues[i])
        end_wingbox = 0.75*c(HalfspanValues[i])
        z_positions = []
        z_diff = []
        y_upper_positions = []
        y_upper_diff = []
        y_lower_positions = []
        y_lower_diff = []
        ds_upper = []
        s_upper = []
        ds_lower = []
        s_lower = []
        
        for i in range(int(N)):
            #print(z)
            if z >= start_wingbox and z <= end_wingbox:
                
                z_positions.append(z)
                y_upper_positions.append(Polyfit_airfoil_upper(z))
                y_lower_positions.append(Polyfit_airfoil_lower(z))
                z = z+dz
             
            else:
                z = z + dz
        
        #print(y_upper_positions)
            
        for i in range(len(z_positions)-1):
            z_diff.append(z_positions[i+1]-z_positions[i])
            y_upper_diff.append(y_upper_positions[i+1]-y_upper_positions[i])
            y_lower_diff.append(y_lower_positions[i+1]-y_lower_positions[i])
        for i in range(len(z_diff)):
            ds_upper.append(m.sqrt(z_diff[i]**2 + y_upper_diff[i]**2))
            ds_lower.append(m.sqrt(z_diff[i]**2 + y_lower_diff[i]**2))
            #print(ds_upper)
        
        s_upper.append(sum(ds_upper))
        s_lower.append(sum(ds_lower))
        #print(s_upper)
        #print(s_lower)
        I_yy = []
        
#        for i in range(len(HalfspanValues)):
#            Iyy = sp.integrate.quad(Polyfit_airfoil_upper(data_z_all_sec[i])*t, 0, s_upper)[0]
#            I_yy.append(Iyy)
#        print(I_yy[0])

#print(airfoil_geometry())