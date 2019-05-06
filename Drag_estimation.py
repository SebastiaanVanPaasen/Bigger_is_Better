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

def Exposed_Areas(root_chord, span, S, taper_ratio, wingtype):
    projected_length_fuselage = D_fuselage/2*(np.tan(LE_sweep)+root_chord*taper_ratio*2/span)
    total_area = projected_length_fuselage*D_fuselage/2
    sweep_excluding_areas = .5*(D_fuselage/2)**2*(np.tan(LE_sweep)+projected_length_fuselage/D_fuselage*2-root_chord)
    exposed_area = S - 2*(total_area - sweep_excluding_areas)
    if wingtype == "Vertical":
        exposed_area= exposed_area/2
    return exposed_area

def Wetted_Areas(D_fuselage, L1, L2, L3):
    S_wing_wet = 1.07*2* Exposed_Areas(root_chord,span,S,taper_ratio, "Horizontal")
    S_horizontal_tail_wet = 1.05 * 2* Exposed_Areas(root_chord_h,span_h,S_h,taper_ratio_h, "Horizontal")
    S_vertical_tail_wet = 1.05 * 2 *Exposed_Areas(root_chord_v,span_v,S_v,taper_ratio_v, "Vertical")
    S_fuselage_wet = np.pi*D_fuselage/4 * (1/(3*L1**2)*((4*L1**2+D_fuselage**2/4)**1.5)\
                                             - D_fuselage + 4*L2 + 2*np.sqrt(L3**2+D_fuselage**2/4))
    wetted_areas = [S_wing_wet, S_horizontal_tail_wet, S_vertical_tail_wet, S_fuselage_wet]
    return (wetted_areas)

def Zero_Lift_Drag_est(S):
    CD_c = [0.0030, 0.0024, 0.006, 0.0025]
    wetted_areas = Wetted_Areas(D_fuselage, L1, L2, L3)
    S_ref = S
    drag_sum = 0
    for i in range(len(CD_c)):
        drag_sum = drag_sum + CD_c[i]*wetted_areas[i]
    CD_0 = (1/S_ref*drag_sum)*1.1
    return CD_0

#print(Zero_Lift_Drag_est(S))

def Drag_est(LE_sweep, A, S, M_cruise, half_chord_sweep, quarter_chord_sweep,alpha,alpha_0,Cl_des):
    e = 4.61*(1-0.045*A**0.68)*np.cos(LE_sweep)**0.15 - 3.1
    K = 1/(np.pi*A*e)
    CD = Zero_Lift_Drag_est(S) + K *\
    Clean_Wing_Lift(A, M_cruise, half_chord_sweep, quarter_chord_sweep, alpha, alpha_0, Cl_des)[1]**2
    return CD 
print("CD is equal to", Drag_est(LE_sweep, A, S, M_cruise, half_chord_sweep, quarter_chord_sweep,alpha,alpha_0,Cl_des))   

CL_CD = Clean_Wing_Lift(A, M_cruise, half_chord_sweep, quarter_chord_sweep, alpha, alpha_0, Cl_des)[1]/Drag_est(LE_sweep, A, S, M_cruise, half_chord_sweep, quarter_chord_sweep,alpha,alpha_0,Cl_des)
print("CL/CD is equal to", CL_CD)