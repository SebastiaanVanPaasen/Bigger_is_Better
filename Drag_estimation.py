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
S_wing_exposed = 500
S_horizontal_tail_exposed = 100
S_vertical_tail_exposed = 50
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

def Wetted_Areas(S_wing_exposed, S_horizontal_tail_exposed, S_vertical_tail_exposed, D_fuselage, L1, L2, L3):
    S_wing_wet = 1.07*2* S_wing_exposed
    S_horizontal_tail_wet = 1.05 * 2* S_horizontal_tail_exposed
    S_vertical_tail_wet = 1.05 * 2 * S_vertical_tail_exposed
    S_fuselage_wet = np.pi*D_fuselage/4 * (1/(3*L1**2)*((4*L1**2+D_fuselage**2/4)**1.5)\
                                             - D_fuselage + 4*L2 + 2*np.sqrt(L3**2+D_fuselage**2/4))
    wetted_areas = [S_wing_wet, S_horizontal_tail_wet, S_vertical_tail_wet, S_fuselage_wet]
    return (wetted_areas)


def Zero_Lift_Drag_est(S):
    CD_c = [0.0030, 0.0024, 0.006, 0.0025]
    wetted_areas = Wetted_Areas(S_wing_exposed, S_horizontal_tail_exposed, S_vertical_tail_exposed, D_fuselage, L1, L2, L3)
    S_ref = S
    drag_sum = 0
    for i in range(len(CD_c)):
        drag_sum = drag_sum + CD_c[i]*wetted_areas[i]
    CD_0 = (1/S_ref*drag_sum)*1.1
    return CD_0

print(Zero_Lift_Drag_est(S))

def Drag_est(LE_sweep, A, S, M_cruise, half_chord_sweep, quarter_chord_sweep,alpha,alpha_0,Cl_des):
    e = 4.61*(1-0.045*A**0.68)*np.cos(LE_sweep)**0.15 - 3.1
    K = 1/(np.pi*A*e)
    CD = Zero_Lift_Drag_est(S) + K *\
    Clean_Wing_Lift(A, M_cruise, half_chord_sweep, quarter_chord_sweep, alpha, alpha_0, Cl_des)[1]**2
    return CD 
print(Drag_est(LE_sweep, A, S, M_cruise, half_chord_sweep, quarter_chord_sweep,alpha,alpha_0,Cl_des))   

CL_CD = Clean_Wing_Lift(A, M_cruise, half_chord_sweep, quarter_chord_sweep, alpha, alpha_0, Cl_des)[1]/Drag_est(LE_sweep, A, S, M_cruise, half_chord_sweep, quarter_chord_sweep,alpha,alpha_0,Cl_des)
print(CL_CD)