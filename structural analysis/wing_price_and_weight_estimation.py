# -*- coding: utf-8 -*-
"""
Created on Fri May 17 10:03:19 2019

@author: mathi
"""
from deflection import s_airfoil
from loading_and_moment_diagrams import load_diagrams
from wingspars import I_zz_spars

def wing_price_weight(N,t_skin,b,length,t_spar):
    density_mat_skin = 2801  # kg/m3
    mat_cost_skin = 6.6  # euro/kg
    density_mat_spar = 2801  # kg/m3
    mat_cost_spar = 6.6
    skin_volume = sum(s_airfoil(N,b))*t_skin
    skin_mass = skin_volume * density_mat_skin
    skin_price = skin_mass * mat_cost_skin
    spar_areas = I_zz_spars(N,b,length,t_spar)[1]+I_zz_spars(N,b,length,t_spar)[2]
    spar_volume = spar_areas*b
    spar_mass = spar_volume * density_mat_spar
    spar_price = spar_mass * mat_cost_spar
    
    total_price = skin_price + spar_price 
    total_mass = skin_mass + spar_mass 
    return skin_mass, spar_mass, total_mass, total_price, 


print(wing_price_weight(100,0.005,62.3,0.3,0.025)[:])


def required_Izz(Cr):
    M_max = max(load_diagrams(100)[0])
    y_max = 0.07 * Cr
    # print(y_max)
    sigma_ult = 552 * 10 ** 6
    I_zz = M_max * y_max / sigma_ult
    print(I_zz)
    return I_zz