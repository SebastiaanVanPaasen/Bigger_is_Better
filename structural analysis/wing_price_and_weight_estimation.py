# -*- coding: utf-8 -*-
"""
Created on Fri May 17 10:03:19 2019

@author: mathi
"""
from deflection import s_airfoil
from loading_and_moment_diagrams import load_diagrams
from Airfoil_inertia import inertia

def wing_price_weight(N,t_skin,b,length,t_spar):
    density_mat_skin = 2801  # kg/m3
    mat_cost_skin = 6.6  # euro/kg
    density_mat_spar = 2801  # kg/m3
    mat_cost_spar = 6.6
    skin_volume = sum(s_airfoil(N,b))*t_skin
    skin_mass = skin_volume * density_mat_skin
    skin_price = skin_mass * mat_cost_skin
    spar_areas = I_zz_spars(N,b,length,t_spar)[1]+I_zz_spars(N,b,length,t_spar)[2]
    spar_volume = spar_areas*b/2
    spar_mass = spar_volume * density_mat_spar
    spar_price = spar_mass * mat_cost_spar
    
    total_price = (skin_price + spar_price)*2
    total_mass = (skin_mass + spar_mass)*2
    return 'skin_mass', skin_mass, 'spar_mass', spar_mass, 'total_mass', total_mass,'total_price', total_price 

t_skin = 0.002
N = 100
t_spar = 0.06
b = 60.#47.83#39.56#41.76
length = 0.3

print(wing_price_weight(N,t_skin,b,length,t_spar)[:])


def I_zz_spars(N,b,length,t_spar):
    
    dy_frontspar = inertia(N,b)[0][0] #upper_y[0] 
    dy_backspar = inertia(N,b)[0][-1] #upper_y[-1]
#    print(dy_frontspar)
    
#    t = 0.025
#    length = 0.3
    area = t_spar*length
#    print('frontsparlength',inertia(N,b)[0][0]-inertia(N,b)[1][0])
#    print('backsparlength',inertia(N,b)[0][-1]-inertia(N,b)[1][-1])

    front_spar_area = t_spar*length*2 + (inertia(N,b)[0][0]-inertia(N,b)[1][0])*t_spar
    back_spar_area = t_spar*length*2 + (inertia(N,b)[0][-1]-inertia(N,b)[1][-1])*t_spar
    
    I_zz_front = ((1/12) * length*t_spar**3 + area*dy_frontspar**2)*2 + (1/12)*t_spar*(inertia(N,b)[0][0]-inertia(N,b)[1][0])**3
    I_zz_back = ((1/12) * length*t_spar**3 + area*dy_backspar**2)*2 + (1/12)*t_spar*(inertia(N,b)[0][-1]-inertia(N,b)[1][-1])**3
    I_zz = I_zz_front + I_zz_back

    return 'I_zz_spars=', I_zz, front_spar_area, back_spar_area   

length = 0.3
t_spar = 0.06
print(I_zz_spars(N,b,length,t_spar))

