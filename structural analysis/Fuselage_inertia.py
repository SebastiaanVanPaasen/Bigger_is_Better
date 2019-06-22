# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 16:10:44 2019

@author: mathi
"""

from thickness_distr_fuselage import t_fus_due_to_moment

x_pos = []
y_pos = []
file = open("C:/Users/mathi/Documents/DSE/Bigger_is_Better/structural analysis/fus_coor.txt", 'r')
for line in file:
    x_pos.append([(float(line.split()[0]))])
    y_pos.append([(float(line.split()[1]))])
x_pos = np.array(x_pos)/1000
y_pos = np.array(y_pos)/1000

def fuselage_centroid(boom_area, y_pos, x_pos):
    
    t_max, alpha_segment, R, segment_length, theta = t_fus_due_to_moment(sigma, Mx, My)


    y_centroid = sum(boom_area*y_pos)/sum(boom_area)
    x_centroid = sum(boom_area*x_pos)/sum(boom_area)
    
    return y_centroid, x_centroid
    
set_boom_area = 0.02

def fuselage_booms(set_boom_area, sigma, Mx, My):
    
    t_max, alpha_segment, R, segment_length, theta = t_fus_due_to_moment(sigma, Mx, My)
    
    req_boom_area = t_max*segment_length
      
    return req_boom_area

    
def fuselage_inertia(boom_area, y_pos, x_pos, My, Mx, theta):
    
    I_xx = sum(boom_area*(y_pos)**2)
    I_yy = sum(boom_area*(x_pos)**2)
    
    sigma_final = ((Mx * np.cos(theta[i]) + My * np.sin(theta[i])) * R_2[j][0] * np.sin(alpha_segment[j]) / I_xx + ( My * np.cos(theta[i]) - Mx * np.sin(theta[i]))*R_2[j][0]*np.cos(alpha_segment[j])/I_yy)/sigma


    return I_xx, I_yy

def fuselage_stiffener_cost(density,cost):
    #option 1: having different boom areas
    req_boom_area = fuselage_booms()
    stiffener_volume = req_boom_area*stiffener_length
    stiffener_mass = stiffener_volume*density
    stiffener_cost = stiffener_mass*cost
    #option 2: having the same boom areas but different pitch
    
    
    return
    