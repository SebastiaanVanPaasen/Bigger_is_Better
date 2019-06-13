# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 16:10:44 2019

@author: mathi
"""

from thickness_distr_fuselage import t_fus_due_to_moment


def fuselage_centroid():
    y_centroid = sum(req_boom_area*R_2*np.sin(alpha_segment))/sum(req_boom_area)
    x_centroid = sum(req_boom_area*x_pos*np.cos(alpha_segment))/sum(req_boom_area)
    
    return y_centroid, x_centroid
    
set_boom_area = 0.02
def fuselage_booms(set_boom_area):
    t_iter, alpha_segment, R_2, segment_length, x_pos = t_fus_due_to_moment()
    req_boom_area = t_iter*segment_length
    n_sec = np.zeros((1,len(req_boom_area)))
    
    for i in range(len(req_boom_area)):
        n_sec = int(req_boom_area[i]/set_boom_area)
#    req_boom_area = req_boom_area - 0.10*segment_length        
        
    return req_boom_area

    
def fuselage_inertia():
    req_boom_area = fuselage_booms()
    I_xx = sum(req_boom_area*(R_2*np.sin(alpha_segment))**2)
    I_yy = sum(req_boom_area*(x_pos*np.cos(alpha_segment))**2)

    return

def fuselage_stiffener_cost(density,cost):
    #option 1: having different boom areas
    req_boom_area = fuselage_booms()
    stiffener_volume = req_boom_area*stiffener_length
    stiffener_mass = stiffener_volume*density
    stiffener_cost = stiffener_mass*cost
    #option 2: having the same boom areas but different pitch
    
    
    return
    