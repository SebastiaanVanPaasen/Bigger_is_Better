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
#qcsweep =   quarter chord sweep
#hcsweeph =   half chord sweep horizontal tail
#hcsweepw = half chord sweep wing
#b_f     =   fuselage width
#b_n     =   nacelle diameter
#S       =   surface area of the wing
#S_net   =   surface area minus the part crossing the fuselage
#h_f     =   fuselage height
#l_fn    =   distance between nose to LE root of the wing
#c       =   chord
#b       =   span
#c_r      =   root chord
#l_h     =   distance between 1/4c wing to horiz. stabilizer
#l_n     =   distance between 1/4c wing to furthest point of nacelle (positive front of 1/4c)
#h_wh     =   vertical distance between horizontal tail and wing
#s_wTEh   =   horizontal distance between horizontal tail and wing

"""Other"""
#M_h_cruise      =   mach based on maximum speed during cruise that horizontal tail faces
#M_c_cruise      =   mach based on maximum speed during cruise that wing faces
#kn              =   reference value 
#V_h_V           =   V_h/V ratio
#eta             =   a constant 
#phi             =   angle between horizontal chord plane and wing chord plane
#n_engine        =   number of engines
#SM              =   stability margin

import numpy as np
from constants_and_conversions import *
from scipy.interpolate import interp1d
import scipy as sp
import matplotlib.pyplot as plt


x_ac_wing = 0.4 #this value needs to be read of a graph which requires (M at max cruise speed, AR, taper, sweep)
SM = 0.05
x_cg = np.linspace(0,1,100)

def C_L_alpha_h(M_h_cruise, eta, hcsweeph, A_h):
    beta = np.sqrt(1-M_h_cruise**2)
    C_L_alpha_h = 2*np.pi*A_h/(2 + np.sqrt(4 + (A_h*beta/eta)*(1+ np.tan(hcsweeph)**2/beta**2)))
    
    return C_L_alpha_h

def C_L_alpha_w(M_w_cruise, eta, hcsweepw, A_w):
    beta = np.sqrt(1-M_w_cruise**2)
    C_L_alpha_w = 2*np.pi*A_w/(2 + np.sqrt(4 + (A_w*beta/eta)*(1+ np.tan(hcsweepw)**2/beta**2)))
    
    return C_L_alpha_w

def C_L_alpha_Ah(M_w_cruise, eta, hcsweepw, A_w, b_f, b, S, c_r):
    S_net = S - (b_f*c_r)
    C_L_alpha_Ah = C_L_alpha_w(M_w_cruise, eta, hcsweepw, A_w) * (1 + 2.15*(b_f/b))*(S_net/S) + (np.pi/2)*(b_f**2/S)

    return C_L_alpha_Ah

def x_ac(b_n_1, b_n_2, b_n_3, b_n_4, l_n_1, l_n_2, l_n_3, l_n_4, x_ac_wing, M_w_cruise, eta, hcsweepw, A_w, b_f, b, c_r, S, h_f, l_fn, c, c_g, qcsweep, taper):
    
    S_net = S - (b_f*c_r)
    x_ac_fus1 = -(1.8/C_L_alpha_Ah(M_w_cruise, eta, hcsweepw, A_w, b_f, b, S, c_r)) *(b_f*h_f*l_fn/(S*c))
    x_ac_fus2 = (0.273/(1+taper))*((b_f*c_g*(b - b_f))/(c*c*(b+2.15*b_f)))*np.tan(qcsweep)
    
    b_n = [b_n_1, b_n_2, b_n_3, b_n_4]
    l_n = [l_n_1, l_n_2, l_n_3, l_n_4]
    x_ac_nac = []
    x_ac_nac_c = 0.
    for i in range(len(b_n)):
        x_ac_nac_c = (b_n[i]*b_n[i]*l_n[i])/(S*c*C_L_alpha_Ah(M_w_cruise, eta, hcsweepw, A_w, b_f, b, S, c_r))
        x_ac_nac.append(x_ac_nac_c)
        
    x_ac = x_ac_wing + x_ac_fus1 + x_ac_fus2 + sum(x_ac_nac)
    
    return x_ac

def de_da(l_h, b, qcsweep, phi, h_wh, s_wTEh, A_w, M_w_cruise, eta, hcsweepw):
    
    r = 2*l_h/b
    w = h_wh**2 + s_wTEh**2
    m_tv = (2/b) * w * np.sin(phi) #(zie plaatje vraag Mathilde)
    k_e_sweep = (0.1124+ 0.1265*qcsweep + 0.1766*qcsweep**2)/r**2 + 0.1024/r +2.
    k_e_sweep_0 = 0.1124/(r*r) + 0.1024/r +2. 
    de_da = (k_e_sweep/k_e_sweep_0)*((r/(r**2+m_tv**2))*(0.4876/np.sqrt(r**2+0.6319 + m_tv**2))+(1+(r**2/(r**2 + 0.7915 + 5.0734*m_tv**2))**(0.3113))*(1 - np.sqrt(m_tv**2/(1+m_tv**2))))*(C_L_alpha_w(M_w_cruise, eta, hcsweepw, A_w))
    
    #de_da_check = 4./(A +2.) #sanity check
    #print("Sanity check = ", de_da_check)
    
    return de_da

def Sh_S_stability(x_cg, M_h_cruise, eta, hcsweeph, hcsweepw, A_w, A_h, M_w_cruise, b_f, b, S, l_h, qcsweep, phi, h_wh, s_wTEh, V_h_V,b_n_1, b_n_2, b_n_3, b_n_4, l_n_1, l_n_2, l_n_3, l_n_4, x_ac_wing, h_f, l_fn, c, c_r, c_g, taper, SM):
    
    S_net = S - (b_f*c_r)
    Sh_S_stability = [] #with SM
    Sh_S_stability_lessSM = []
    for i in range(len(x_cg)):
        den = (C_L_alpha_h(M_h_cruise, eta, hcsweeph, A_h)/C_L_alpha_Ah(M_w_cruise, eta, hcsweepw, A_w, b_f, b, S, c_r))*(1-de_da(l_h, b, qcsweep, phi, h_wh, s_wTEh, A_w, M_w_cruise, eta, hcsweepw))*(l_h/c)*(V_h_V)
        Sh_S = (1/den)*x_cg[i] - (x_ac(b_n_1, b_n_2, b_n_3, b_n_4, l_n_1, l_n_2, l_n_3, l_n_4, x_ac_wing, M_w_cruise, eta, hcsweepw, A_w, b_f, b, c_r, S, h_f, l_fn, c, c_g, qcsweep, taper) - SM)/den
        Sh_S_stability.append(Sh_S)
        den = (C_L_alpha_h(M_h_cruise, eta, hcsweeph, A_h)/C_L_alpha_Ah(M_w_cruise, eta, hcsweepw, A_w, b_f, b, S, c_r))*(1-de_da(l_h, b, qcsweep, phi, h_wh, s_wTEh, A_w, M_w_cruise, eta, hcsweepw))*(l_h/c)*(V_h_V)
        Sh_S_less = (1/den)*x_cg[i] - x_ac(b_n_1, b_n_2, b_n_3, b_n_4, l_n_1, l_n_2, l_n_3, l_n_4, x_ac_wing, M_w_cruise, eta, hcsweepw, A_w, b_f, b, c_r, S, h_f, l_fn, c, c_g, qcsweep, taper)/den
        Sh_S_stability_lessSM.append(Sh_S_less)
    
    stability_trend = np.polyfit(x_cg, Sh_S_stability, 1)
#    print(stability_trend[0])
#    print(stability_trend[1])
#    plt.plot(x_cg, Sh_S_stability, 'r')
#    plt.plot(x_cg, Sh_S_stability_lessSM)
#    plt.ylim(bottom=0)
#    plt.show()
#    
    return Sh_S_stability

