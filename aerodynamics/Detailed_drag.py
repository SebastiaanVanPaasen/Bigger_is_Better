#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 12:40:59 2019

@author: Max
"""
import numpy as np
#Detailed drag esstiation


rho = 1.225
V = 200
mu = 1
M = 0.7
k = 0.634*10**(-5)


#rho is cruise density
#V is cruise speed
#l is the characteristic lenth
#mu is the viscostity for cruise conditions
#k surface finish parameter smooth paint should be selected
#M is the mach number
#trans is the transition point  (ratio between 0 and 1)


def friction_coefficient(rho,V, l, mu, k, M, trans ):
    Re = min(((rho*V*l)/(mu)),38.21*(l/k)**(1.053))
    Cf_lam = (1.328)/(np.sqrt(Re))
    Cf_tur = (0.455)/((np.log10(Re))**(2.58)*(1+0.144*M**2)**(0.65))
    Cf = trans*Cf_lam + (1-trans)*Cf_tur
    return Cf

#form factor calculations  
#FF for lifting bodys like wing, tail, strut and pylon
def FF_lifting(max_t_loc,t_c, M,sweep_m ):
    FF = (1 + (0.6/max_t_loc)*t_c + 100*t_c**100)*(1.34*M**0.18*(np.cos(sweep_m))**0.28)
    return FF

#FF for fuselage and cannopy
def FF_body(l,d):
    f = l/d
    FF = 1 + 60/(f**3) + f/400
    return FF

#FF for nacelles and external storage
def FF_nacelle(l,d):
    f = l/d
    FF = 1 + 0.35/f
    return FF

#calulate wetted areas
    
def Fus_wetted_area(d, L1, L2, L3):
    wetted_area = ((np.pi*d)/4)*((1/(3*L1**2))*(((4*L1**2 + (d**2)/4)**1.15) - (d**3)/8) - d + 4*L2 + 2*np.sqrt(L3**2 + (d**2)/4))
    return wetted_area


#Set interference factors for each component see ADSEE-II lecture 2 slide 48
IF_wing = 1.1
IF_emp = 1.05
IF_fus = 1.0
IF_nacelle = 1.1
IF_strut = 1.1
    

    