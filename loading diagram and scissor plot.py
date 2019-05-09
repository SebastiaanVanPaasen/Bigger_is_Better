# -*- coding: utf-8 -*-
"""
Created on Thu May  9 08:38:50 2019

@author: Hidde
"""

import numpy as np
import matplotlib.pyplot as plt
from class_I_empennage_landinggear import *
from fuselage_cross_section import *


# ---------------INPUTS--------------------------------------------------------

g              = 9.81    # [m/s^2]
W_payload      = 500000  # [N]
n_cargo        = 2       # Number of cargo compartments (1 or 2)
m_passenger    = 90      # mass of passenger and hand baggage kg

cargo_fwdfrac  = 0.6     # Fraction of the amount the front compartment holds, if n_cargo = 1 then the value is 1

# Weights from class II in [N]
W_fuse      = 80000  
W_nlg       = 1500
W_vt        = 2500
W_ht        = 3500
W_fix       = 50000
W_wing      = 70000
W_nac       = 8000
W_prop      = 15000
W_mlg       = 6000

# Import cg values from class I sizing

xcg_fuse, xcg_emp, xcg_fix, xcg_nac, xcg_prop, xcg_wing, S_h_frac_S, S_v_frac_S, X_LEMAC, xcg_fwd, xcg_aft, zcg = classI_empennage(MAC, l_fuselage, x_eng, l_n, xcg_OEW_MAC, mass_frac_OEW, xcg_payload, mass_frac_payload, xcg_fuel, mass_frac_fuel, D_fuse, b)

xcg_nlg = l_nosecone  # assumption, use landing gear function !!!
xcg_mlg = 0.75*MAC    # assumption, use landing gear function !!!

# Inputs from fuselage dimensions



#------------------------------------------------------------------------------

def potato(W_fuse, W_nlg, W_vt, W_ht, W_fix, W_wing, W_nac, W_prop, W_mlg, xcg_fuse, xcg_emp, xcg_fix, xcg_nac, xcg_prop, xcg_wing, xcg_nlg, xcg_mlg, X_LEMAC, W_payload, n_cargo, l_fuselage, MAC, n_pax, m_passenger, g, cargo_fwdfrac):

    # Fuselage group weight and cg
    W_fusegroup   = W_fuse + W_nlg + W_vt + W_ht + W_fix
    xcg_vt        = xcg_emp
    xcg_ht        = xcg_emp
    xcg_fusegroup = (xcg_fuse*W_fuse + xcg_nlg*W_nlg + xcg_vt*W_vt + xcg_ht*W_ht + xcg_fix*W_fix)/W_fusegroup
    
    # Wing group weight and cg
    W_winggroup   = W_wing + W_nac + W_prop + W_mlg
    xcg_winggroup = (xcg_wing*W_wing + xcg_nac*W_nac + xcg_prop*W_prop + xcg_mlg*W_mlg)/W_winggroup
    xcg_wg_MAC    = xcg_winggroup/MAC
    
    # Aircraft OEW weight and cg
    xcg_OEW_MAC = (xcg_fusegroup + W_winggroup/W_fusegroup*xcg_wg_MAC - X_LEMAC)/(1. + W_winggroup/W_fusegroup) 
    W_OEW       = W_fusegroup + W_winggroup
    
    #------------Cargo loading-------------------------------------------------
    W_pax = n_pax*m_passenger*g
    W_cargo = W_payload - W_pax
    
    if n_cargo == 1:
        xcg_cargo1_MAC  = (0.45*l_fuselage - X_LEMAC)/MAC   # cg cargo as percentage of MAC
        xcg_cargo2_MAC  = 0
        W_cargo1        = W_cargo         
        W_cargo2        = 0
        
    if n_cargo == 2:
        xcg_cargo1_MAC  = (0.3*l_fuselage - X_LEMAC)/MAC   # cg cargo 1 as percentage of MAC
        xcg_cargo2_MAC  = (0.6*l_fuselage - X_LEMAC)/MAC   # cg cargo 2 as percentage of MAC
        W_cargo1        = W_cargo*(cargo_fwdfrac)          
        W_cargo2        = W_cargo*(1-cargo_fwdfrac)        
    
    #--------------- Backward loading------------------------------------------
    
    # Weight order
    W_Orderback      = [W_OEW] + [W_cargo2, W_cargo1]
    
    # cg order
    xcg_Orderback    = [xcg_OEW_MAC] + [xcg_cargo2_MAC, xcg_cargo1_MAC] 
    
    # Weight actual values
    W_Backloading     = [W_OEW]
    for i in range(1, len(W_Orderback)):
        W_Backloading.append(W_Backloading[i-1]+W_Orderback[i])
    
    # cg actual values
    xcg_Backloading   = [xcg_OEW_MAC]
    for i in range(1, len(xcg_Orderback)):
        xcg_Backloading.append((xcg_Backloading[i-1]*W_Backloading[i-1]+xcg_Orderback[i]*W_Orderback[i])/W_Backloading[i])
    
    #----------------Forward loading-------------------------------------------
    
    # Weight order
    W_Orderfwd      = [W_OEW] + [W_cargo1, W_cargo2]
    
    # cg order
    xcg_Orderfwd    = [xcg_OEW_MAC] + [xcg_cargo1_MAC, xcg_cargo2_MAC]
    
    # Weight actual values
    W_Fwdloading     = [W_OEW]
    for i in range(1, len(W_Orderfwd)):
        W_Fwdloading.append(W_Fwdloading[i-1]+W_Orderfwd[i])
    
    # cg actual values
    xcg_Fwdloading   = [xcg_OEW_MAC]
    for i in range(1, len(xcg_Orderfwd)):
        xcg_Fwdloading.append((xcg_Fwdloading[i-1]*W_Fwdloading[i-1]+xcg_Orderfwd[i]*W_Orderfwd[i])/W_Fwdloading[i])
    
    
    print(W_OEW, W_cargo1, W_cargo2, xcg_OEW_MAC, xcg_cargo1_MAC, xcg_cargo2_MAC, W_Fwdloading, xcg_Fwdloading)
    
    plt.plot(xcg_Backloading, W_Backloading)
    plt.plot(xcg_Fwdloading, W_Fwdloading)
    plt.show()
    print(X_LEMAC)
    return(xcg_fusegroup, xcg_winggroup, xcg_OEW_MAC)
    

    
potato(W_fuse, W_nlg, W_vt, W_ht, W_fix, W_wing, W_nac, W_prop, W_mlg, xcg_fuse, xcg_emp, xcg_fix, xcg_nac, xcg_prop, xcg_wing, xcg_nlg, xcg_mlg, X_LEMAC, W_payload, n_cargo, l_fuselage, MAC, n_pax, m_passenger, g, cargo_fwdfrac)
