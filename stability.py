# -*- coding: utf-8 -*-
"""
Created on Thu May  9 14:08:37 2019

@author: Mathilde
"""
"""Geometry"""
#
#A_w     =   aspect ratio of the wing
#A_h     =   aspect ratio of the horizontal tail
#taper   =   taper ratio of the wing
#qcsweep =   quarter chord sweep [rad]
#hcsweep =   half chord sweep [rad]
#b_f     =   fuselage width
#b_n     =   nacelle diameter
#S       =   surface area of the wing
#S_net   =   surface area minus the part crossing the fuselage
#h_f     =   fuselage height
#l_fn    =   distance between nose to LE root of the wing
#c       =   chord
#b       =   span
#l_h     =   distance between 1/4c wing to horiz. stabilizer
#l_n     =   distance between 1/4c wing to furthest point of nacelle (positive front of 1/4c)
#
"""Other"""
#M_h_cruise      =   mach based on maximum speed during cruise that horizontal tail faces
#M_c_cruise      =   mach based on maximum speed during cruise that wing faces
#k_n              =   reference value 
#V_h_V           =   V_h/V ratio
#eta             =   a constant 
#phi             =   angle between horizontal chord and wing chord 
#n_engine        =   number of engines


import numpy as np

x_ac_wing = 0.4 #this value needs to be read of a graph which requires (M at max cruise speed, AR, taper, sweep)

def C_L_alpha_h(M_h_cruise, eta, hcsweep, A_h):
    beta = np.sqrt(1-M_h_cruise**2)
    C_L_alpha_h = 2*pi*A_h/(2 + np.sqrt(4 + (A_h*beta/eta)*(1+ tan(hcsweep)**2/beta**2)))
    
    return C_L_alpha_h

def C_L_alpha_w(M_w_cruise, eta, hcsweep, A_h):
    beta = np.sqrt(1-M_w_cruise**2)
    C_L_alpha_w = 2*pi*A_w/(2 + np.sqrt(4 + (A_w*beta/eta)*(1+ tan(hcsweep)**2/beta**2)))
    
    return C_L_alpha_w

def C_L_alpha_Ah(M_w_cruise, eta, hcsweep, A_h, b_f, b, S_net, S):
    
    C_L_alpha_Ah = C_L_alpha_w(M_w_cruise, eta, hcsweep, A_h) * (1 + 2.15*(b_f/b))*(S_net/S) + (pi/2)*(b_f**2/S)

    return C_L_alpha_Ah

def x_ac(b_n_1, b_n_2, b_n_3, b_n_4, l_n_1, l_n_2, l_n_3, l_n_4, x_ac_wing, M_w_cruise, eta, hcsweep, A_h, b_f, b, S_net, S, h_f, l_fn, S, c, c_g, taper, qcsweep):
    x_ac_fus1 = -(1.8/C_L_alpha_Ah(M_w_cruise, eta, hcsweep, A_h, b_f, b, S_net, S)) *(b_f*h_f*l_fn/(S*c))
    x_ac_fus2 = (0.273/(1+taper))*((b_f*c_g*(b - b_f))/(c*c*(b+2.15*b_f)))*tan(qcsweep)
    b_n = [b_n_1, b_n_2, b_n_3, b_n_4]
    l_n = [l_n_1, l_n_2, l_n_3, l_n_4]
    x_ac_nac = []
    x_ac_nac_c = 0.
    for i in range(len(b_n)):
        x_ac_nac_c = (b_n*b_n*l_n)/(S*c*C_L_alpha_Ah(M_w_cruise, eta, hcsweep, A_h, b_f, b, S_net, S))
        x_ac_nac.append(x_ac_nac_c)
        
    x_ac = x_ac_wing + x_ac_fus1 + x_ac_fus2 + sum(x_ac_nac)
    
    return x_ac

def de_da():
    k_e_sweep = 0.1124 + 0.1265*qcsweep +0.1766*qcsweep**2
    
    

