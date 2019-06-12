# -*- coding: utf-8 -*-
"""
Created on Wed May 15 16:01:47 2019

@author: floyd
"""
import numpy as np
import math as m
import scipy as sp
from loading_and_moment_diagrams import c
from parameter_requirements import *
#import loading_and_moment_diagrams as lm
from airfoil_geometry import *
import centroid_wing as cw
from spar_locations import spar_loc
#from stress_distribution_wing import wing_stress

N = 100 
b = lm.b#60.#47.83#39.56#41.76
HalfspanValues = np.linspace(0, b / 2 - 0.00001, N)
data_z_all_sec = airfoil_geometry(N,b)[0]
data_y_upper_all_sec = airfoil_geometry(N,b)[1]
data_y_lower_all_sec = airfoil_geometry(N,b)[2]
z_centroid_all_sec, y_centroid_all_sec, y_loc_spar_up, y_loc_spar_low, y_loc_stiff_up, y_loc_stiff_low, y_vertical_spar, z_loc_stiff_up, spar_loc_sec, z_loc_stiff_low = wing_centroid(cw.boom_area, cw.spar_areas_hori, cw.t_spar_v, cw.z_c_airfoil, cw.y_c_airfoil, cw.n_stiff_up, cw.n_stiff_low, cw.HalfspanValues)

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
    
    
def I_zz_spars(l_spar_h ,t_spar_v, t_spar_h):

    I_zz_spar = np.zeros((len(HalfspanValues), 1))
    I_yy_spar = np.zeros((len(HalfspanValues), 1))
    I_yz_spar = np.zeros((len(HalfspanValues), 1))
    
    for j in range(len(HalfspanValues)):
        for i in range(len(spar_areas_hori)):
            I_zz_up = (1/12) * l_spar_h*t_spar_h**3 + l_spar_h*t_spar_h*(y_loc_spar_up[j][i]-y_centroid_all_sec[j])**2  
            I_yy_up = (1/12) * t_spar_h*l_spar_h**3 + l_spar_h*t_spar_h*(spar_loc_sec[j][i]-z_centroid_all_sec[j])**2
            I_yz_up = l_spar_h*t_spar_h*(y_loc_spar_up[j][i]-y_centroid_all_sec[j])*(spar_loc_sec[j][i]-z_centroid_all_sec[j])
            
            I_zz_low = (1/12) * l_spar_h*t_spar_h**3 + l_spar_h*t_spar_h*(y_centroid_all_sec[j] - y_loc_spar_low[j][i])**2
            I_yy_low = (1/12) * t_spar_h*l_spar_h**3 + l_spar_h*t_spar_h*(spar_loc_sec[j][i]-z_centroid_all_sec[j])**2
            I_yz_low = l_spar_h*t_spar_h*(y_loc_spar_low[j][i]-y_centroid_all_sec[j])*(spar_loc_sec[j][i]-z_centroid_all_sec[j])
            
            l_spar_v = y_loc_spar_up[j][i] - y_loc_spar_low[j][i]
#            print("l_spar_v", l_spar_v)
            I_zz_vertical = (1/12) * t_spar_v* l_spar_v**3 + l_spar_v*t_spar_v*(y_vertical_spar[j][i] - y_centroid_all_sec[j])**2
            I_yy_vertical = (1/12) * l_spar_v* t_spar_v**3 + l_spar_v*t_spar_v*(spar_loc_sec[j][i] - z_centroid_all_sec[j])**2
            I_yz_vertical = l_spar_v* t_spar_v*(y_vertical_spar[j][i] - y_centroid_all_sec[j])*(spar_loc_sec[j][i] - z_centroid_all_sec[j])
            
            I_zz = I_zz_up + I_zz_low + I_zz_vertical
            I_yy = I_yy_up + I_yy_low + I_yy_vertical
            I_yz = I_yz_up + I_yz_low + I_yz_vertical
            
            I_zz_spar[j] = I_zz
            I_yy_spar[j] = I_yy
            I_yz_spar[j] = I_yz

    return I_zz_spar, I_yy_spar, I_yz_spar

#print(I_zz_spars(cw.l_spar_h, cw.t_spar_v, cw.t_spar_h)[0][0])
#print(I_zz_spars(cw.l_spar_h, cw.t_spar_v, cw.t_spar_h)[1][0])
#print(I_zz_spars(cw.l_spar_h, cw.t_spar_v, cw.t_spar_h)[2][0])

I_zz_spars = I_zz_spars(cw.l_spar_h, cw.t_spar_v, cw.t_spar_h)
I_zz_req = required_Izz(HalfspanValues, data_y_all_sec, y_centroid_all_sec, Mz)

#first define the wing lay out that is required to obtain the right moment of inertia
def wing_geometry(I_zz_req, I_zz_spars, y_loc_stiff_up, y_loc_stiff_low):
    
    sum_boom_areas = np.zeros((len(HalfspanValues), 1))
    single_boom_area = np.zeros((len(HalfspanValues), 1))
    
    for i in range(len(HalfspanValues)):
        
        y_2 = 0
        
        for k in range(len(y_loc_stiff_up[0])):
            y_2 += y_loc_stiff_up[i][k]**2 
    
        for j in range(len(y_loc_stiff_low[0])):
            y_2 += y_loc_stiff_low[i][j]**2 
        
        sum_boom_areas[i] = (I_zz_req[i][0] - I_zz_spars[i][0])/y_2
    
        single_boom_area[i] = sum_boom_areas[i]/(cw.n_stiff_low + cw.n_stiff_up)
    
    plt.plot(HalfspanValues, single_boom_area)
    plt.show()
    return sum_boom_areas, single_boom_area

#print("single boom area",wing_geometry(I_zz_req, I_zz_spars, y_loc_stiff_up, y_loc_stiff_low)[1]) 
I_zz_spar = I_zz_spars[0]
I_yy_spar = I_zz_spars[1]
I_yz_spar = I_zz_spars[2]
airfoil_area = cw.airfoil_area
z_c_airfoil = cw.z_c_airfoil
y_c_airfoil = cw.y_c_airfoil
#calculate the wing moment of inertia per section 
def inertia_wing(I_zz_spar, I_yy_spar, I_yz_spar, airfoil_area, z_c_airfoil, y_c_airfoil):
    
    I_zz = np.zeros((len(HalfspanValues), 1))
    I_yy = np.zeros((len(HalfspanValues), 1))
    I_yz = np.zeros((len(HalfspanValues), 1))
    
    
    for i in range(len(HalfspanValues)):
        
        I_zz_booms=0
        I_yy_booms=0
        I_yz_booms=0
        I_zz_airfoil = 0
        I_yy_airfoil = 0
        I_yz_airfoil = 0
        
        for j in range(len(y_loc_stiff_up[0])):
            
            I_zz_booms += boom_area*(y_loc_stiff_up[i][j] - y_centroid_all_sec[i])**2
            I_yy_booms += boom_area*(z_loc_stiff_up[i][j] - z_centroid_all_sec[i])**2
            I_yz_booms += boom_area*(y_loc_stiff_up[i][j] - y_centroid_all_sec[i])*(z_loc_stiff_up[i][j] - z_centroid_all_sec[i])
        
        for k in range(len(y_loc_stiff_low[0])):
            
            I_zz_booms += boom_area*(y_centroid_all_sec[i]- y_loc_stiff_low[i][k])**2
            I_yy_booms += boom_area*(z_loc_stiff_low[i][k] - z_centroid_all_sec[i])**2
            I_yz_booms += boom_area*(y_centroid_all_sec[i]-y_loc_stiff_low[i][k])*(z_loc_stiff_low[i][k] - z_centroid_all_sec[i])
        
        I_zz_airfoil += airfoil_area[i]*(y_c_airfoil[i] - y_centroid_all_sec[i])**2
        I_yy_airfoil += airfoil_area[i]*(z_c_airfoil[i] - z_centroid_all_sec[i])**2
        I_yz_airfoil += airfoil_area[i]*(y_c_airfoil[i] - y_centroid_all_sec[i])*(z_c_airfoil[i] - z_centroid_all_sec[i])
        
        I_zz[i] = I_zz_booms + I_zz_spar[i][0] + I_zz_airfoil
        I_yy[i] = I_yy_booms + I_yy_spar[i][0] + I_yy_airfoil
        I_yz[i] = I_yz_booms + I_yz_spar[i][0] + I_yz_airfoil
    
    return I_zz, I_yy, I_yz
    
print(inertia_wing(I_zz_spar, I_yy_spar, I_yz_spar, airfoil_area, z_c_airfoil, y_c_airfoil)[0])
print(inertia_wing(I_zz_spar, I_yy_spar, I_yz_spar, airfoil_area, z_c_airfoil, y_c_airfoil)[1][0])
print(inertia_wing(I_zz_spar, I_yy_spar, I_yz_spar, airfoil_area, z_c_airfoil, y_c_airfoil)[2][0])
    
#    defining the integrated function

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
        
        
#    upper_y = np.zeros(len(spar_loc))
#    lower_y = np.zeros(len(spar_loc))
#    
#    for i in range(len(spar_loc)):
#        
#        upper_y[i] = Polyfit_airfoil_upper(spar_loc_sec[x_position][i])
#    
#        lower_y[i] = Polyfit_airfoil_lower(spar_loc_sec[x_position][i])
#        
#    delta_z = []
#    delta_y_upper = []
#    delta_y_lower = []
#    middle_y_upper = []
#    middle_y_lower = []
#    
#    for i in range(len(spar_loc)-1):
#        delta_y_upper.append(upper_y[i+1]-upper_y[i])
#        delta_y_lower.append(lower_y[i+1] - lower_y[i])
#        delta_z.append(spar_loc_sec[x_position][i+1]-spar_loc_sec[x_position][i])
#        middle_y_upper.append(upper_y[i] + (upper_y[i+1]-upper_y[i])/2)
#        middle_y_lower.append(lower_y[i] + (lower_y[i+1]-lower_y[i])/2)
#    
#    middle_z = spar_loc_sec[x_position] + delta_z[0]
#    
#    #Moment of inertia of flanges on its own centroid
#    I_xx_upper = []    
#    I_xx_lower = []
#    I_yy_upper = []
#    I_yy_lower = []
#    I_xy_upper = []
#    I_xy_lower = []
#    
#    upper_s = []
#    lower_s = []
#    
#    t= 0.005
#    centroid_y = 0.05
#    centroid_z = 0.4*c(HalfspanValues[0])
#    
#    for i in range(len(delta_z)):
#        s_upper = np.sqrt(delta_z[i]**2 + delta_y_upper[i]**2)
#        s_lower = np.sqrt(delta_z[i]**2 + delta_y_lower[i]**2)
#        upper_s.append(s_upper)
#        lower_s.append(s_lower)
#        beta = m.atan(delta_y_upper[i]/delta_z[i])
#        area_upper = s_upper*t
#        area_lower = s_lower*t
#        dy_upper = middle_y_upper[i] - centroid_y
#        dy_lower = middle_y_lower[i] - centroid_y
#        dx_lower = middle_z[i] - centroid_z
#        dx_upper = middle_z[i] - centroid_z
#        I_xx_lower.append(s_lower**3*t*np.sin(beta)**2/12.+ area_lower*dy_lower**2)
#        I_xx_upper.append(s_upper**3*t*np.sin(beta)**2/12.+ area_upper*dy_upper**2)
#        I_yy_lower.append(s_lower**3*t*np.cos(beta)**2/12.+ area_lower*dx_lower**2)
#        I_yy_upper.append(s_upper**3*t*np.cos(beta)**2/12. + area_upper*dx_upper**2)
#        I_xy_lower.append(s_lower**3*t*np.sin(2*beta)/24.+ area_lower*dx_lower*dy_lower)
#        I_xy_upper.append(s_upper**3*t*np.sin(2*beta)/24. + area_upper*dx_upper*dy_upper)
#    
#    return upper_y, lower_y, upper_s, lower_s#y_centroid_sec, z_centroid_sec

#y_centroid_sec, z_centroid_sec = inertia(N,b)
#upper_y, lower_y, upper_s, lower_s
#print(sum(inertia(N,b)[2]))
#sum_s_wingbox = sum(inertia(N,b)[2]) + sum(inertia(N,b)[3])
#print(sum_s_wingbox)
#print(sum(I_xx_lower)+sum(I_xx_lower) + sum(I_xx_upper) + sum(I_xx_upper))


