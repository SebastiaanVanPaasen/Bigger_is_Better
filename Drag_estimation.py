# -*- coding: utf-8 -*-
"""
Created on Mon May  6 15:01:51 2019

@author: floyd
"""
import numpy as np
from Lift_estimation import Clean_Wing_Lift

D_fuselage  = 3
L1 = 5
L2 = 30
L3 = 5
S = 500
A = 10
M_cruise = 0.7
half_chord_sweep = 0.5
quarter_chord_sweep = 0.6
LE_sweep = 0.65
alpha = 0.5
alpha_0 = 0.3
Cl_des = 1.2
airfoil_eff_factor = 0.95
root_chord = 15
span = np.sqrt(S*A)
taper_ratio = 0.3
root_chord_h = 5
root_chord_v = 4
span_h = 20
span_v = 25
taper_ratio_h = 0.5
taper_ratio_v = 0.5
S_h = 50
S_v = 70
tip_chord = root_chord*taper_ratio

#def Exposed_Areas(root_chord, span, S, taper_ratio, wingtype): #exposed wing area calulation
#    projected_length_fuselage = D_fuselage/2*(np.tan(LE_sweep)+root_chord*taper_ratio*2/span)
#    total_area = projected_length_fuselage*D_fuselage/2
#    sweep_excluding_areas = .5*(D_fuselage/2)**2*(np.tan(LE_sweep)+projected_length_fuselage/D_fuselage*2-root_chord)
#    exposed_area = S - 2*(total_area - sweep_excluding_areas)
#    if wingtype == "Vertical":
#        exposed_area= exposed_area/2
#    return exposed_area

def Wing_wetted_area(root_chord, tip_chord, D_fuselage, span, S):
    fuselage_chord = root_chord*(1-(D_fuselage/span)) #determine the cord at the outer fuselage
    exposed_area = S-(0.5*(fuselage_chord+root_chord))*D_fuselage #calulate the exposed area, substracting the area inside the plane from the S
    wetted_area = 1.07*2*exposed_area #wetted area aproximation
    return wetted_area
    
def H_tail_wetted_area(root_chord_h, taper_ratio_h, span_h):
    exposed_area = 0.5*(root_chord_h*(1+taper_ratio_h))*span_h
    wetted_area = 1.05*2*exposed_area
    return wetted_area

def V_tail_wetted_area(root_chord_v, taper_ratio_v, span_v):
    exposed_area = 0.5*(root_chord_v*(1+taper_ratio_v))*(span_v/2)
    wetted_area = 1.05*2*exposed_area
    return wetted_area

def Fus_wetted_area(D_fuselage, L1, L2, L3):
    wetted_area = np.pi*D_fuselage/4 * (1/(3*L1**2)*((4*L1**2+D_fuselage**2/4)**1.5)\
                                           - D_fuselage + 4*L2 + 2*np.sqrt(L3**2+D_fuselage**2/4))
    return wetted_area

def Zero_Lift_Drag_est(S, Wing_wetted_area, H_wetted_area, V_wetted_area, Fus_wetted_area):
    CD_0 = (1/S)*(0.003*Wing_wetted_area+ 0.0024*Fus_wetted_area + 0.0025*(H_wetted_area+V_wetted_area))
    return CD_0

Wing_wetted_area = Wing_wetted_area(root_chord, tip_chord, D_fuselage, span, S)
H_wetted_area = H_tail_wetted_area(root_chord_h, taper_ratio_h, span_h)
V_wetted_area = V_tail_wetted_area(root_chord_v, taper_ratio_v, span_v)
Fus_wetted_area = Fus_wetted_area(D_fuselage, L1, L2, L3)
CD_0 = Zero_Lift_Drag_est(S, Wing_wetted_area, H_wetted_area, V_wetted_area, Fus_wetted_area)

print(Wing_wetted_area,H_wetted_area, V_wetted_area, Fus_wetted_area, CD_0)

