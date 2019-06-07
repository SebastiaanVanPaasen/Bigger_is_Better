# -*- coding: utf-8 -*-
"""
Created on Wed May 15 16:01:47 2019

@author: floyd
"""
import numpy as np
import math as m
import scipy as sp
from stress_distribution_wing import load_airfoil
from loading_and_moment_diagrams import c
import loading_and_moment_diagrams as lm
from airfoil_geometry import airfoil_geometry

N = 100 
b = lm.b#60.#47.83#39.56#41.76
HalfspanValues = np.linspace(0, b / 2 - 0.00001, N)


def s_airfoil(N,b):
    
    data_z_all_sec = airfoil_geometry(N,b)[0]
    data_y_upper_all_sec = airfoil_geometry(N,b)[1]
    data_y_lower_all_sec = airfoil_geometry(N,b)[2]
    ds_sec_all = []
    s_all_sec = []
    
    for i in range(len(data_z_all_sec)):
        for j in range(len(data_z_all_sec[0])-1):
            ds_sec = np.sqrt((data_z_all_sec[i][j+1]-data_z_all_sec[i][j])**2+(data_y_upper_all_sec[i][j+1]-data_y_upper_all_sec[i][j])**2) + np.sqrt((data_z_all_sec[i][j+1]-data_z_all_sec[i][j])**2+(data_y_lower_all_sec[i][j+1]-data_y_lower_all_sec[i][j])**2)
            ds_sec_all.append(ds_sec)
    
    ds_sec_all = np.asarray(ds_sec_all)
    ds_sec_all = np.reshape(ds_sec_all, (len(data_z_all_sec),len(data_z_all_sec[0])-1))
    
    for i in range(len(ds_sec_all)):
        s_all_sec.append(sum(ds_sec_all[i]))
  
    return(s_all_sec)
    
    


def inertia(N,b):
    #loading the airfoil data
#    data = load_airfoil("naca3414.txt")
#    dataz = np.asarray(data[1])
#    datay = np.asarray(data[2])
#    #ordering the data
    x_position = 0
    data_z_all_sec = airfoil_geometry(N,b)[0]
    data_y_upper_all_sec = airfoil_geometry(N,b)[1]
    data_y_lower_all_sec = airfoil_geometry(N,b)[2]
    
#    defining the integrated function
    Polyfit_airfoil_upper = sp.interpolate.interp1d(data_z_all_sec[x_position], data_y_upper_all_sec[x_position], kind="cubic", fill_value="extrapolate")
    Polyfit_airfoil_lower = sp.interpolate.interp1d(data_z_all_sec[x_position], data_y_lower_all_sec[x_position], kind="cubic", fill_value="extrapolate")

    spar_loc_sec = []
    for i in range(len(HalfspanValues)):
        nr_spars = 4
        first_spar_location = 0.2*c(HalfspanValues[i])
        last_spar_location = 0.75*c(HalfspanValues[i])
        delta_spar = (last_spar_location-first_spar_location)/(nr_spars)
        spar_loc = np.arange(first_spar_location,last_spar_location+delta_spar,delta_spar)
        spar_loc_sec.append(spar_loc)
#
#    print(spar_loc_sec)        
    
    #centroid berekenen:
#    y_centroid_sec = []
#    z_centroid_sec = []
#    t = 0.01
#    for i in range(1,len(data_y_upper_all_sec)):
#        y_centroid_point = []
#        z_centroid_point = []
#        ds_sec_all = 0
#        for j in range(len(data_y_upper_all_sec[0])):
#            y_cen_inertia_above =((abs(data_y_upper_all_sec[i][j]-data_y_upper_all_sec[i][j-1])/2.)+data_y_upper_all_sec[i][j-1])
#            y_cen_inertia_below =((abs(data_y_lower_all_sec[i][j]-data_y_lower_all_sec[i][j-1])/2.)+data_y_lower_all_sec[i][j-1])
#            z_cen_inertia = ((abs(data_z_all_sec[i][j]-data_z_all_sec[i][j-1])/2.)+data_z_all_sec[i][j-1])
#            print(y_cen_inertia_below)
#            ds_airfoil_upper_sec = np.sqrt(abs(data_z_all_sec[i][j]-data_z_all_sec[i][j-1])**2+ abs(data_y_upper_all_sec[i][j]-data_y_upper_all_sec[i][j-1])**2)
#            ds_airfoil_lower_sec = np.sqrt(abs(data_z_all_sec[i][j]-data_z_all_sec[i][j-1])**2+ abs(data_y_lower_all_sec[i][j]-data_y_lower_all_sec[i][j-1])**2)
#            y_centroid_above_sec = ds_airfoil_upper_sec*t*y_cen_inertia_above
#            y_centroid_below_sec = ds_airfoil_lower_sec*t*y_cen_inertia_below
#            z_centroid_above_sec = ds_airfoil_upper_sec*t*z_cen_inertia
#            z_centroid_below_sec = ds_airfoil_lower_sec*t*z_cen_inertia
#            y_centroid_point.append(y_centroid_above_sec)
#            y_centroid_point.append(y_centroid_below_sec)
#            z_centroid_point.append(z_centroid_above_sec)
#            z_centroid_point.append(z_centroid_below_sec)
#            ds_sec_all = ds_sec_all+ds_airfoil_upper_sec+ds_airfoil_lower_sec
#        y_centroid = sum(y_centroid_point)/(ds_sec_all*t)
#        z_centroid = sum(z_centroid_point)/(ds_sec_all*t)
#        y_centroid_sec.append(y_centroid)
#        z_centroid_sec.append(z_centroid)
    #print(y_centroid)
    #print(z_centroid)
#    plt.plot(data_z_all_sec[-1],data_y_upper_all_sec[-1])
#    plt.plot(data_z_all_sec[-1],data_y_lower_all_sec[-1])
#    plt.scatter(z_centroid_sec[-1],y_centroid_sec[-1])
#    plt.show()
        
        
    upper_y = np.zeros(len(spar_loc))
    lower_y = np.zeros(len(spar_loc))
    
    for i in range(len(spar_loc)):
        
        upper_y[i] = Polyfit_airfoil_upper(spar_loc_sec[x_position][i])
    
        lower_y[i] = Polyfit_airfoil_lower(spar_loc_sec[x_position][i])
        
    delta_z = []
    delta_y_upper = []
    delta_y_lower = []
    middle_y_upper = []
    middle_y_lower = []
    
    for i in range(len(spar_loc)-1):
        delta_y_upper.append(upper_y[i+1]-upper_y[i])
        delta_y_lower.append(lower_y[i+1] - lower_y[i])
        delta_z.append(spar_loc_sec[x_position][i+1]-spar_loc_sec[x_position][i])
        middle_y_upper.append(upper_y[i] + (upper_y[i+1]-upper_y[i])/2)
        middle_y_lower.append(lower_y[i] + (lower_y[i+1]-lower_y[i])/2)
    
    middle_z = spar_loc_sec[x_position] + delta_z[0]
    
    #Moment of inertia of flanges on its own centroid
    I_xx_upper = []    
    I_xx_lower = []
    I_yy_upper = []
    I_yy_lower = []
    I_xy_upper = []
    I_xy_lower = []
    
    upper_s = []
    lower_s = []
    
    t= 0.005
    centroid_y = 0.05
    centroid_z = 0.4*c(HalfspanValues[0])
    
    for i in range(len(delta_z)):
        s_upper = np.sqrt(delta_z[i]**2 + delta_y_upper[i]**2)
        s_lower = np.sqrt(delta_z[i]**2 + delta_y_lower[i]**2)
        upper_s.append(s_upper)
        lower_s.append(s_lower)
        beta = m.atan(delta_y_upper[i]/delta_z[i])
        area_upper = s_upper*t
        area_lower = s_lower*t
        dy_upper = middle_y_upper[i] - centroid_y
        dy_lower = middle_y_lower[i] - centroid_y
        dx_lower = middle_z[i] - centroid_z
        dx_upper = middle_z[i] - centroid_z
        I_xx_lower.append(s_lower**3*t*np.sin(beta)**2/12.+ area_lower*dy_lower**2)
        I_xx_upper.append(s_upper**3*t*np.sin(beta)**2/12.+ area_upper*dy_upper**2)
        I_yy_lower.append(s_lower**3*t*np.cos(beta)**2/12.+ area_lower*dx_lower**2)
        I_yy_upper.append(s_upper**3*t*np.cos(beta)**2/12. + area_upper*dx_upper**2)
        I_xy_lower.append(s_lower**3*t*np.sin(2*beta)/24.+ area_lower*dx_lower*dy_lower)
        I_xy_upper.append(s_upper**3*t*np.sin(2*beta)/24. + area_upper*dx_upper*dy_upper)
    
    return upper_y, lower_y, upper_s, lower_s#y_centroid_sec, z_centroid_sec

#y_centroid_sec, z_centroid_sec = inertia(N,b)
#upper_y, lower_y, upper_s, lower_s
#print(sum(inertia(N,b)[2]))
#sum_s_wingbox = sum(inertia(N,b)[2]) + sum(inertia(N,b)[3])
#print(sum_s_wingbox)
#print(sum(I_xx_lower)+sum(I_xx_lower) + sum(I_xx_upper) + sum(I_xx_upper))


