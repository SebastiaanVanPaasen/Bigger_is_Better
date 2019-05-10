# -*- coding: utf-8 -*-
"""
Created on Thu May  9 08:38:50 2019

@author: Hidde
"""

import numpy as np
import matplotlib.pyplot as plt
from class_I_empennage_landinggear import *
from fuselage_cross_section import *
from seats import *


# ---------------INPUTS--------------------------------------------------------

g              = 9.81    # [m/s^2]
W_payload      = 450000 # [N]
n_cargo        = 2       # Number of cargo compartments (1 or 2)
m_passenger    = 90      # mass of passenger and hand baggage kg

cargo_fwdfrac  = 0.6     # Fraction of the amount the front compartment holds, if n_cargo = 1 then the value is 1

# Weights from class II in [N]
W_fuse      = 80000*5  
W_nlg       = 1500*5
W_vt        = 2500*5
W_ht        = 3500*5
W_fix       = 50000*5
W_wing      = 70000*5
W_nac       = 8000*5
W_prop      = 15000*5
W_mlg       = 6000*5

W_fuel      = 938780

xcg_fuel    = 0.55   # % MAC

Safety_margin = 0.02   # for cg range

Npax = 330
Npax_below = 330
d_inner, d_outer, lcabin, lcabin_below, lcabin_above, l_tailcone_range, l_nosecone_range, l_fuselage,tot_seating_abreast, N_aisle_above, N_aisle_below,N_rows_above,N_rows_below,lpax_below,lpax_above,tot_seating_abreast_last_row,Pseat=fuselage_cross_section(Npax,Npax_below)
xcg_seats = cg_seats(Pseat, tot_seating_abreast,N_rows_above,N_rows_below,lpax_below,lpax_above,l_nosecone_range,tot_seating_abreast_last_row)
W_window, W_aisle, W_middle = W_seats(tot_seating_abreast,N_rows_above,N_rows_below,lpax_below,lpax_above,l_nosecone_range,tot_seating_abreast_last_row)
D_fuse = d_outer
l_nosecone = (l_nosecone_range[0] + l_nosecone_range[1]) / 2.
xcg_payload = l_fuselage * 0.5  # cg-location payload w.r.t. nose [m]
# Import cg values from class I sizing

cg_locations, S_h_frac_S, S_v_frac_S, X_LEMAC, x_h = class_I_empennage(MAC, l_fuselage, x_eng, l_n, xcg_OEW_MAC, mass_frac_OEW, xcg_payload, mass_frac_payload, xcg_fuel, mass_frac_fuel, D_fuse, b)

xcg_nlg = l_nosecone  # assumption, use landing gear function !!!
xcg_mlg = 0.75*MAC    # assumption, use landing gear function !!!

# Inputs from fuselage dimensions





#------------------------------------------------------------------------------

def potato(Safety_margin, xcg_fuel, W_fuel, W_fuse, W_nlg, W_vt, W_ht, W_fix, W_wing, W_nac, W_prop, W_mlg, cg_locations, xcg_nlg, xcg_mlg, X_LEMAC, W_payload, n_cargo, l_fuselage, MAC, n_pax, m_passenger, g, cargo_fwdfrac, xcg_seat, W_window, W_aisle, W_middle):
    
    xcg_fuse = cg_locations[0]
    xcg_emp = cg_locations[1]
    xcg_fix = cg_locations[2]
    xcg_nac = cg_locations[3]
    xcg_prop = cg_locations[4]
    xcg_wing = cg_locations[5]
    # Fuselage group weight and cg
    #W_fusegroup   = W_fuse + W_nlg + W_vt + W_ht + W_fix
    W_fusegroup = 2866210*0.4*0.35
    xcg_vt        = xcg_emp
    xcg_ht        = xcg_emp
    xcg_fusegroup = (xcg_fuse*W_fuse + xcg_nlg*W_nlg + xcg_vt*W_vt + xcg_ht*W_ht + xcg_fix*W_fix)/W_fusegroup
    
    # Wing group weight and cg
    #W_winggroup   = W_wing + W_nac + W_prop + W_mlg
    W_winggroup = 2866210*0.4*0.17
    xcg_winggroup = (xcg_wing*W_wing + xcg_nac*W_nac + xcg_prop*W_prop + xcg_mlg*W_mlg)/W_winggroup
    xcg_wg_MAC    = xcg_winggroup/MAC
    
    min_cg = []
    max_cg = []
    X_LEMAC_range = []

    for j in np.arange((X_LEMAC/l_fuselage)-0.3, (X_LEMAC/l_fuselage)+0.3, 0.001):
        
        X_LEMAC_range.append(j)
        X_LEMAC = j*l_fuselage
        
        # Aircraft OEW weight and cg
        #xcg_OEW_MAC = (xcg_fusegroup + (W_winggroup/W_fusegroup)*xcg_wg_MAC - X_LEMAC)/(1. + W_winggroup/W_fusegroup) 
        W_OEW       = W_fusegroup + W_winggroup
        xcg_OEW_MAC = ((xcg_fusegroup-X_LEMAC)/MAC*W_fusegroup + xcg_wg_MAC*W_winggroup)/W_OEW
        
        
        
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
        W_Orderback      = [W_OEW] + [W_cargo2, W_cargo1]+list((np.array(W_window[::-1])*m_passenger*g))+list((np.array(W_aisle[::-1])*m_passenger*g))+list((np.array(W_middle[::-1])*m_passenger*g)) + [W_fuel]
        
        # cg order
        xcg_Orderback    = [xcg_OEW_MAC] + [xcg_cargo2_MAC, xcg_cargo1_MAC]+3*list((np.array(xcg_seats[::-1])-X_LEMAC)/MAC) + [xcg_fuel]
        
        
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
        W_Orderfwd      = [W_OEW] + [W_cargo1, W_cargo2] + list((np.array(W_window)*m_passenger*g))+list((np.array(W_aisle)*m_passenger*g))+list((np.array(W_middle)*m_passenger*g)) + [W_fuel]
        
        # cg order
        xcg_Orderfwd    = [xcg_OEW_MAC] + [xcg_cargo1_MAC, xcg_cargo2_MAC] + 3*list((np.array(xcg_seats)-X_LEMAC)/MAC) + [xcg_fuel]
        
        # Weight actual values
        W_Fwdloading     = [W_OEW]
        for i in range(1, len(W_Orderfwd)):
            W_Fwdloading.append(W_Fwdloading[i-1]+W_Orderfwd[i])
        
        # cg actual values
        xcg_Fwdloading   = [xcg_OEW_MAC]
        for i in range(1, len(xcg_Orderfwd)):
            xcg_Fwdloading.append((xcg_Fwdloading[i-1]*W_Fwdloading[i-1]+xcg_Orderfwd[i]*W_Orderfwd[i])/W_Fwdloading[i])
        
        
        min_cg.append(min(xcg_Fwdloading)-Safety_margin)
        max_cg.append(max(xcg_Backloading)+Safety_margin)
      
#        plt.plot(xcg_Backloading, W_Backloading)
#        plt.plot(xcg_Fwdloading, W_Fwdloading)
#        plt.show()
#       
#    plt.plot(min_cg, X_LEMAC_range)
#    plt.plot(max_cg, X_LEMAC_range)
#    plt.show()
    return(xcg_fusegroup, xcg_winggroup, xcg_OEW_MAC,min_cg, max_cg,X_LEMAC_range)
    

    
potato(Safety_margin, xcg_fuel, W_fuel, W_fuse, W_nlg, W_vt, W_ht, W_fix, W_wing, W_nac, W_prop, W_mlg, cg_locations, xcg_nlg, xcg_mlg, X_LEMAC, W_payload, n_cargo, l_fuselage, MAC, n_pax, m_passenger, g, cargo_fwdfrac, xcg_seats, W_window, W_aisle, W_middle)
