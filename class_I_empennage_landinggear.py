# -*- coding: utf-8 -*-
"""
Created on Mon May  6 13:32:31 2019

@author: Hidde
"""

from fuselage_cross_section import *


# Inputs

# MAC = 8.  # mean aerodynamic chord [m]
# n_pax = 450  # Number of passengers
# n_paxbelow = 450
# x_eng = -1.  # x-location engines w.r.t. XLEMAC [m]
# l_n = 2.  # length nacelle [m]
# xcg_OEW_MAC = 0.25  # initial cg location OEW w.r.t. MAC [m]
# mass_frac_OEW = 0.5  # mass fraction OEW
# mass_frac_payload = 0.2  # mass fraction payload
# xcg_fuel = 22  # cg-location fuel w.r.t nose [m]
# mass_frac_fuel = 0.3  # mass fraction fuel
# b = 30.  # wingspan [m]
# MLW = 1500000  # maximum landing weight [N]
# MTOW = 2000000  # maximum take-off weight [N]
#
# d_inner, d_outer, lcabin, lcabin_below, lcabin_above, l_tailcone_range, l_nosecone_range, l_fuselage, tot_seating_abreast, N_aisle_above, N_aisle_below = fuselage_cross_section(
#     n_pax, n_paxbelow)
#
# D_fuse = d_outer
# l_nosecone = (l_nosecone_range[0] + l_nosecone_range[1]) / 2.
# xcg_payload = l_fuselage * 0.5  # cg-location payload w.r.t. nose [m]

mass_frac_wing = 0.117
mass_frac_emp = 0.023
mass_frac_fuse = 0.098
mass_frac_nac = 0.018
mass_frac_prop = 0.072
mass_frac_fix = 0.118
m_fracs = [mass_frac_wing, mass_frac_emp, mass_frac_fuse, mass_frac_nac, mass_frac_prop, mass_frac_fix, 0.3, 0.3, 0.4]

def class_I_empennage(MAC, l_fuselage, x_eng, l_n, xcg_OEW_MAC, mass_frac_OEW, xcg_payload, mass_frac_payload, xcg_fuel,
                      mass_frac_fuel, D_fuse, b):
    # Mass fractions main components, landing gear neglected
    mass_frac_wing = 0.117
    mass_frac_emp = 0.023
    mass_frac_fuse = 0.098
    mass_frac_nac = 0.018
    mass_frac_prop = 0.072
    mass_frac_fix = 0.118

    # xcg-locations for wing mounted engines aircraft

    # Wing group - in [m] w.r.t. X-LEMAC
    xcg_wing = 0.4 * MAC
    xcg_prop = x_eng + 0.4 * l_n
    xcg_nac = x_eng + 0.4 * l_n

    # Fuselage group - in [m] w.r.t. nose
    xcg_fuse = 0.4 * l_fuselage
    xcg_fix = 0.4 * l_fuselage
    xcg_emp = 0.9 * l_fuselage

    # Averaged mass fractions and cg locations

    # Wing group
    mass_frac_winggroup = mass_frac_wing + mass_frac_prop + mass_frac_nac
    xcg_winggroup = (
                            xcg_wing * mass_frac_wing + xcg_prop * mass_frac_prop + xcg_nac * mass_frac_nac) / mass_frac_winggroup

    # Fuselage group
    mass_frac_fusegroup = mass_frac_fuse + mass_frac_fix + mass_frac_emp
    xcg_fusegroup = (
                            xcg_fuse * mass_frac_fuse + xcg_fix * mass_frac_fix + xcg_emp * mass_frac_emp) / mass_frac_fusegroup

    # X location of the leading edge of the mean aerodynamic chord
    X_LEMAC = xcg_fusegroup + MAC * (
            (xcg_winggroup / MAC) * (mass_frac_winggroup / mass_frac_fusegroup) - xcg_OEW_MAC * (
            1 + (mass_frac_winggroup / mass_frac_fusegroup)))

    # cg excursion
    xcg_OEW = xcg_OEW_MAC * MAC + X_LEMAC
    xcglist = []
    xcglist.append(xcg_OEW)
    xcglist.append((xcg_OEW * mass_frac_OEW + xcg_payload * mass_frac_payload) / (mass_frac_OEW + mass_frac_payload))
    xcglist.append((xcg_OEW * mass_frac_OEW + xcg_fuel * mass_frac_fuel) / (mass_frac_OEW + mass_frac_fuel))
    xcglist.append((xcg_OEW * mass_frac_OEW + xcg_payload * mass_frac_payload + xcg_fuel * mass_frac_fuel) / (
            mass_frac_OEW + mass_frac_payload + mass_frac_fuel))

    # Maximum and minimum cg locations w.r.t. nose [m]
    xcg_fwd = min(xcglist)
    xcg_aft = max(xcglist)

    # Z location cg
    zcg = 0.5 * D_fuse

    # Empennage sizing

    # Assumption: tail arms are taken from the cg of the tail
    # Horizontal and vertical tail locations

    
    

    # Normalized tail coefficients obtained from statistics
    
    
    
    
    S = 300.
    b = 50.
    
    x_h = xcg_emp
    V_h_norm = 0.99
    tap_h            = 0.4
    A_h              = 4.
    sweep_quartchord_h = np.pi*30./180.
    n = 1. 
    print(x_h)
    while n > 0.001:
        
        S_frac_h = (V_h_norm * MAC) / (x_h - xcg_aft)
        S_h = S_frac_h*S
        
        b_h = np.sqrt(A*S_h)
        C_r_h = (2.*S_h)/(b_h*(1.+tap_h))
        C_t_h = C_r_h*tap_h
        MAC_h = (2./3.)*C_r_h*((1.+tap_h+tap_h**2)/(1.+tap_h))
        sweep_LE_h = np.arctan(np.tan(sweep_quartchord_h)-((C_r_h)/(2.*b_h))*(tap_h-1))
        sweep_TE_h = np.arctan((C_t_h-C_r_h)/(b_h/2.)+np.tan(sweep_LE_h))
        y_MAC_h = (b_h/6.)*((1.+2.*tap_h)/(1.+tap_h))
        fuse_MAC_h = MAC_h-0.25*MAC_h-y_MAC_h*np.tan(sweep_TE_h)
        
        n = abs((x_h-(l_fuselage-fuse_MAC_h))/(x_h))
        x_h = l_fuselage-fuse_MAC_h
        
    print(S_h, b_h, C_r_h, C_t_h, MAC_h, fuse_MAC_h, y_MAC_h)
    print(x_h)
    print()
    x_v = xcg_emp    
    V_v_norm = 0.081
    tap_v = 0.4
    A_v = 1.5
    sweep_LE_v = np.pi*35./180.
    n = 1.
    print(x_v)
    while n > 0.001:
        
        S_frac_v = (V_v_norm * b) / (x_v - xcg_aft)
        S_v = S_frac_v*S
        
        b_v = np.sqrt(A*S_v)
        C_r_v = (S_v)/(b_v*(1.+tap_v))
        C_t_v = C_r_v*tap_v
        MAC_v = (2./3.)*C_r_v*((1.+tap_v+tap_v**2)/(1.+tap_v))
        sweep_TE_v = np.arctan((C_t_v-C_r_v)/(b_v)+np.tan(sweep_LE_v))
        y_MAC_v = (b_v/3.)*((1.+2.*tap_v)/(1.+tap_v))
        fuse_MAC_v = MAC_v-0.25*MAC_v-y_MAC_v*np.tan(sweep_TE_v)
        
        n = abs((x_v-(l_fuselage-fuse_MAC_v))/(x_v))
        x_v = l_fuselage-fuse_MAC_v
    print(sweep_TE_v, S_v, b_v, C_r_v, C_t_v, MAC_v, fuse_MAC_v,y_MAC_v)
    print(x_v)
    print()
    print(l_fuselage)
    
    cg_locations = np.array([xcg_fuse, xcg_emp, xcg_fix, xcg_nac, xcg_prop, xcg_wing, xcg_fwd, xcg_aft, zcg])

    return cg_locations, S_h_frac_S, S_v_frac_S, X_LEMAC, x_h, x_v

class_I_empennage(m_fracs, MAC, l_fuselage, x_eng, l_n, xcg_OEW_MAC, xcg_payload, xcg_fuel, D_fuse, b)

# class_I_empennage(MAC, l_fuselage, x_eng, l_n, xcg_OEW_MAC, mass_frac_OEW, xcg_payload, mass_frac_payload, xcg_fuel,
#                  mass_frac_fuel, D_fuse, b)


def size_tail(wing_area, volume_fraction, tail_sweep, aspect_ratio):
    tail_area = volume_fraction * wing_area
    taper_ratio = 0.2 * (2 - tail_sweep * np.pi / 180)
    tail_span = np.sqrt(tail_area * aspect_ratio)
    root_chord = (2 * tail_area) / ((1 + taper_ratio) * tail_span)
    tip_chord = root_chord * taper_ratio
    return root_chord, tip_chord, tail_span, tail_area


def classI_landinggear(MTOW, MLW, l_nosecone):
    N_mw = MLW / 210000
    N_mw = int(N_mw)

    while N_mw % 4 != 0:
        N_mw = N_mw + 1

    if N_mw <= 12:
        N_s = 2
    else:
        N_s = 4

    # l_nlg =

    return N_mw, N_s

# print (classI_landinggear(MTOW, MLW, l_nosecone))
