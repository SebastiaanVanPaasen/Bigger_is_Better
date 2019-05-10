# -*- coding: utf-8 -*-
"""
Created on Mon May  6 13:32:31 2019

@author: Hidde
"""

from fuselage_cross_section import *

# Inputs

MAC = 8.  # mean aerodynamic chord [m]
n_pax = 450  # Number of passengers
n_paxbelow = 450
x_eng = -1.  # x-location engines w.r.t. XLEMAC [m]
l_n = 2.  # length nacelle [m]
xcg_OEW_MAC = 0.25  # initial cg location OEW w.r.t. MAC [m]
mass_frac_OEW = 0.5  # mass fraction OEW
mass_frac_payload = 0.2  # mass fraction payload
xcg_fuel = 22  # cg-location fuel w.r.t nose [m]
mass_frac_fuel = 0.3  # mass fraction fuel
b = 30.  # wingspan [m]
MLW = 1500000  # maximum landing weight [N]
MTOW = 2000000  # maximum take-off weight [N]


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

    x_h = xcg_emp
    x_v = xcg_emp

    # Normalized tail coefficients obtained from statistics
    V_h_norm = 0.99
    V_v_norm = 0.081

    S_h_frac_S = (V_h_norm * MAC) / (x_h - xcg_aft)
    S_v_frac_S = (V_v_norm * b) / (x_v - xcg_aft)
    cg_locations = np.array([xcg_fuse, xcg_emp, xcg_fix, xcg_nac, xcg_prop, xcg_wing, xcg_fwd, xcg_aft, zcg])

    return cg_locations, S_h_frac_S, S_v_frac_S, X_LEMAC, x_h, x_v


def size_tail(wing_area, volume_fraction, tail_sweep, aspect_ratio):
    tail_area = volume_fraction*wing_area
    taper_ratio = 0.2*(2 - tail_sweep*np.pi/180)
    tail_span = np.sqrt(tail_area*aspect_ratio)
    root_chord = (2*tail_area)/((1+taper_ratio)*tail_span)
    tip_chord = root_chord*taper_ratio
    return root_chord, tip_chord, tail_span, tail_area


def classI_landinggear(MTOW, MLW, l_nosecone):
    N_mw = MLW / 210000.
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
