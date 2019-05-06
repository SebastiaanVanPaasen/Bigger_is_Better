#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 12:16:04 2019

@author: Max
"""
import numpy as np #delete later should be loaded in the general program


def planform(M_cruise, CL_cruise, surface_area, aspect_ratio):
    
    #Sweep calulation
    M_t = 0.935 #Technoligy factor for super critical airfoild
    M_dd = M_cruise + 0.03 #drag divergence mach number
    
    #Torenbeek estimation for quarter cord sweep values are in radians
    if M_cruise<=0.7:
        quarter_cord_sweep = 0.
    else:
        quarter_cord_sweep = np.arccos(0.75*(M_t/M_dd))
    #calulate taper ratio using Torenbeek
    taper_ratio = 0.2*(2.-quarter_cord_sweep)
    
    #calulate span, root and tip cords geometrically
    span = np.sqrt(surface_area*aspect_ratio)
    cord_root = (2*surface_area)/((1+taper_ratio)*span)
    cord_tip = taper_ratio*cord_root
    
    #Select dihedral using Adsee aproach
    dihedral = np.deg2rad(3 -(np.rad2deg(quarter_cord_sweep)/10)+2) #selected for low wing if hight wing use -2 instead of +2
    
    #tickness over cord ratio
    leading_edge_sweep = np.arctan(np.tan(quarter_cord_sweep)-(cord_root/(2*span)*(taper_ratio-1)))
    print(leading_edge_sweep)
    half_cord_sweep = np.arctan(((cord_tip/2+span/2*np.tan(leading_edge_sweep))-(cord_root/2))/(span/2))
    tickness_over_cord = (np.cos(half_cord_sweep)**3*(M_t-M_dd*np.cos(half_cord_sweep))-0.115*CL_cruise**1.5)/(np.cos(half_cord_sweep)**2)
    tickness_over_cord = min(tickness_over_cord,0.18)
    
    return(quarter_cord_sweep, taper_ratio, span, cord_root, cord_tip, dihedral, tickness_over_cord)

M_cruise = 0.8      #inputs from different part 
surface_area = 500  #inputs from different part 
aspect_ratio = 9    #inputs from different part 
CL_cruise = 0.5

quarter_cord_sweep, taper_ratio, span, cord_root, cord_tip, dihedral, tickness_over_cord = planform(M_cruise, CL_cruise, surface_area, aspect_ratio)
print(quarter_cord_sweep, taper_ratio, span, cord_root, cord_tip, dihedral, tickness_over_cord)
    
    
    








                       
                       
    

   



