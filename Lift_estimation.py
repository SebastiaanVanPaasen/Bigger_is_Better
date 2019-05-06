# -*- coding: utf-8 -*-
"""
Created on Mon May  6 12:20:48 2019

@author: floyd
"""
import numpy as np
from matplotlib import pyplot as plt
## Lift Estimation
#%%  #VALUES ARE USED FOR FIRST ESTIMATE USE THE ACTUAL INPUTS!
A = 9
M_cruise = 0.7
half_chord_sweep = 0.5
quarter_chord_sweep = 0.6
alpha = 0.5
alpha_0 = 0.3
Cl_des = 1.2
airfoil_eff_factor = 0.95
#%%

# Function to find CL_alpha, CL and alpha trim in the case of no wing twist
def Clean_Wing_Lift(A, M_cruise, half_chord_sweep, quarter_chord_sweep, alpha, alpha_0, Cl_des):
    beta = np.sqrt(1-M_cruise) #prandtl_glauert correction
    #DATCOM method for CL_alpha
    CL_alpha = (2*np.pi* A)/( 2+ np.sqrt(4 + (A * beta/airfoil_eff_factor)**2 *(1+np.tan(half_chord_sweep)**2 / beta**2)))
    CL = CL_alpha*(alpha - alpha_0)                 #calculate the Lift
    CL_des = Cl_des/np.cos(quarter_chord_sweep)**2 #CL_des based on Cl_des from airfoil
    alpha_0_L = alpha_0 #In case of no wing twist
    alpha_trim = CL_des/CL_alpha + alpha_0_L    #Find trim angle to obtain required lift
    return(CL_alpha, CL, alpha_trim)

