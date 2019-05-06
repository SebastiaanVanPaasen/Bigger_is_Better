# -*- coding: utf-8 -*-
"""
Created on Mon May  6 13:32:31 2019

@author: Hidde
"""
# Inputs

MAC               = 4.         # mean aerodynamic chord [m]
l_fuse            = 40.        # length fuselage [m]
x_eng             = -1.        # x-location engines w.r.t. XLEMAC [m]
l_n               = 2.         # length nacelle [m]          
xcg_OEW_MAC       = 0.25       # initial cg location OEW w.r.t. MAC [m]
mass_frac_OEW     = 0.5        # mass fraction OEW
xcg_payload       = 18         # cg-location payload w.r.t. nose [m]  
mass_frac_payload = 0.2        # mass fraction payload
xcg_fuel          = 22         # cg-location fuel w.r.t nose [m]
mass_frac_fuel    = 0.3        # mass fraction fuel 
D_fuse            = 2.         # diameter fuselage [m]
b                 = 30.        # wingspan [m]

def classI_empennage(MAC, l_fuse, x_eng, l_n, xcg_OEW_MAC, mass_frac_OEW, xcg_payload, mass_frac_payload, xcg_fuel, mass_frac_fuel, D_fuse, b):
    # Mass fractions main components, landing gear neglected 
    mass_frac_wing = 0.117
    mass_frac_emp  = 0.023
    mass_frac_fuse = 0.098
    mass_frac_nac  = 0.018
    mass_frac_prop = 0.072
    mass_frac_fix  = 0.118
    
    # xcg-locations for wing mounted engines aircraft
    
    # Wing group - in [m] w.r.t. X-LEMAC
    xcg_wing = 0.4*MAC
    xcg_prop = x_eng + 0.4*l_n
    xcg_nac  = x_eng + 0.4*l_n
    
    # Fuselage group - in [m] w.r.t. nose
    xcg_fuse = 0.4*l_fuse
    xcg_fix  = 0.4*l_fuse
    xcg_emp  = 0.9*l_fuse
    
    # Averaged mass fractions and cg locations
    
    # Wing group
    mass_frac_winggroup = mass_frac_wing + mass_frac_prop + mass_frac_nac
    xcg_winggroup       = (xcg_wing*mass_frac_wing + xcg_prop*mass_frac_prop + xcg_nac*mass_frac_nac)/mass_frac_winggroup
    
    # Fuselage group
    mass_frac_fusegroup = mass_frac_fuse + mass_frac_fix + mass_frac_emp
    xcg_fusegroup       = (xcg_fuse*mass_frac_fuse + xcg_fix*mass_frac_fix + xcg_emp*mass_frac_emp)/mass_frac_fusegroup
    
    # X location of the leading edge of the mean aerodynamic chord
    X_LEMAC = xcg_fusegroup + MAC*((xcg_winggroup/MAC)*(mass_frac_winggroup/mass_frac_fusegroup) - xcg_OEW_MAC*(1 + (mass_frac_winggroup/mass_frac_fusegroup)))
    
    # cg excursion
    xcg_OEW = xcg_OEW_MAC*MAC + X_LEMAC
    xcglist = []
    xcglist.append(xcg_OEW)
    xcglist.append((xcg_OEW*mass_frac_OEW + xcg_payload*mass_frac_payload)/(mass_frac_OEW + mass_frac_payload))
    xcglist.append((xcg_OEW*mass_frac_OEW + xcg_fuel*mass_frac_fuel)/(mass_frac_OEW + mass_frac_fuel))
    xcglist.append((xcg_OEW*mass_frac_OEW + xcg_payload*mass_frac_payload + xcg_fuel*mass_frac_fuel)/(mass_frac_OEW + mass_frac_payload + mass_frac_fuel))
    
    # Maximum and minimum cg locations w.r.t. nose [m]
    xcg_fwd = min(xcglist)
    xcg_aft = max(xcglist)
    
    # Z location cg
    zcg = 0.5*D_fuse
    
    # Empennage sizing
    
    # Assumption: tail arms are taken from the cg of the tail
    # Horizontal and vertical tail locations
    
    x_h = xcg_emp
    x_v = xcg_emp
    
    # Normalized tail coefficients obtained from statistics
    V_h_norm = 0.99
    V_v_norm = 0.081
    
    S_h_frac_S = (V_h_norm*MAC)/(x_h - xcg_aft)
    S_v_frac_S = (V_v_norm*b)/(x_v - xcg_aft)
    
    return(S_h_frac_S, S_v_frac_S, X_LEMAC, xcg_fwd, xcg_aft, zcg)
    
classI_empennage(MAC, l_fuse, x_eng, l_n, xcg_OEW_MAC, mass_frac_OEW, xcg_payload, mass_frac_payload, xcg_fuel, mass_frac_fuel, D_fuse, b)

def landinggearlocation():
    
    return


