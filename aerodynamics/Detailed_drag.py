#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 12:40:59 2019

@author: Max
"""
import numpy as np
#Detailed drag esstiation


rho = 0.466348
a = 303.793
mu = 0.0000150498
M = 0.775
k = 0.634*10**(-5)
V = M*a

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

#FF_strut = FF_lifting(0.3,0.24,M,0)
#Cf_strut = friction_coefficient(rho,V,0.094, mu,k,M,0.1)
#Sw_strut = 4
#S_ref = 183
#
#Cd_strut = (1/S_ref)*FF_strut*Cf_strut*Sw_strut
#print(Cd_strut) 

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
    S_wet = ((np.pi*d)/4)*((1/(3*L1**2))*(((4*L1**2 + (d**2)/4)**1.15) - (d**3)/8) - d + 4*L2 + 2*np.sqrt(L3**2 + (d**2)/4))
    return S_wet

def Wing_wetted_area(c_r, c_tip, d_fus, b, S,):  
    taper = c_t/c_r
    c_f = c_r * (1 - (d_fus / b)*(1-taper))  # determine the cord at the outer fuselage
    exposed_area = S - (0.5 *(c_f + c_r))*d_fus  # calulate the exposed area, substracting the area inside the plane from the S
    print (exposed_area)
    S_wet = 1.07*2*exposed_area  # wetted area approximation
    return S_wet

def Tail_wetted_area(S_exp):
    S_wet = 1.05*2*S_exp
    return S_wet




#S_wet for the nacelle should be determined later source of tu delft available ask max
#Set interference factors for each component see ADSEE-II lecture 2 slide 48
IF_wing = 1.1
IF_emp = 1.05
IF_fus = 1.0
IF_nacelle = 1.1
IF_strut = 1.1

l_fus = 53.286
d_fus = 5.951

L1 = 11.009
L3 = 19.638
L2 = l_fus- L1 - L3

S = 249.281

#l_fus = 51.029
#d_fus = 6.288
#
#L1 = 11.632
#L3 = 20.749
#L2 = l_fus- L1 - L3
#
#S = 247.386

Cf_fus = friction_coefficient(rho,V, l_fus , mu, k, M, 0.0 )
FF_fus = FF_body(l_fus, d_fus)
S_wet_fus = Fus_wetted_area(d_fus, L1, L2, L3)

Cd0_fus = (Cf_fus*FF_fus*S_wet_fus)/S
print (Cf_fus, FF_fus, S_wet_fus)
print (Cd0_fus)
    