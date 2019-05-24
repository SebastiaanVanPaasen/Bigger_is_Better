# -*- coding: utf-8 -*-
"""
Created on Mon May  6 13:32:31 2019

@author: Hidde
"""

from class_I.fuselage_cross_section import *


# Inputs
# used for the iterator
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


def class_I_empennage(mass_frac, MAC, l_fuselage, x_eng, l_n, xcg_OEW_MAC, xcg_payload, xcg_fuel, D_fuse, b, S, taper,
                      v_tail, LE_sweep, h_tail):
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
    mass_frac_winggroup = mass_frac[0] + mass_frac[4] + mass_frac[3]
    xcg_winggroup = (
                            xcg_wing * mass_frac[0] + xcg_prop * mass_frac[4] + xcg_nac * mass_frac[
                        3]) / mass_frac_winggroup

    # Fuselage group
    mass_frac_fusegroup = mass_frac[2] + mass_frac[5] + mass_frac[1]
    xcg_fusegroup = (
                            xcg_fuse * mass_frac[2] + xcg_fix * mass_frac[5] + xcg_emp * mass_frac[
                        1]) / mass_frac_fusegroup

    # X location of the leading edge of the mean aerodynamic chord
    X_LEMAC = xcg_fusegroup + MAC * (
            (xcg_winggroup / MAC) * (mass_frac_winggroup / mass_frac_fusegroup) - xcg_OEW_MAC * (
            1 + (mass_frac_winggroup / mass_frac_fusegroup)))
    Y_MAC = (b / 6.) * ((1 + 2 * taper) / (1 + taper))

    # cg excursion
    xcglist = []
    xcg_OEW = xcg_OEW_MAC * MAC + X_LEMAC
    xcglist.append(xcg_OEW)
    xcglist.append((xcg_OEW * mass_frac[6] + xcg_payload * mass_frac[7]) / (mass_frac[6] + mass_frac[7]))
    xcglist.append((xcg_OEW * mass_frac[6] + xcg_fuel * mass_frac[8]) / (mass_frac[6] + mass_frac[8]))
    xcglist.append((xcg_OEW * mass_frac[6] + xcg_payload * mass_frac[7] + xcg_fuel * mass_frac[8]) / (
            mass_frac[6] + mass_frac[7] + mass_frac[8]))

    # Maximum and minimum cg locations w.r.t. nose [m]
    xcg_fwd = min(xcglist)
    xcg_aft = max(xcglist)

    # Z location cg
    zcg = 0.5 * D_fuse

    # Empennage sizing
    # Assumption: tail arms are taken from the cg of the tail
    # Horizontal and vertical tail locations
    V_h_norm = 0.99
    V_v_norm = 0.081
    A_h, taper_h, QC_sweep_h = h_tail[0], h_tail[1], h_tail[2]
    A_v, taper_v, QC_sweep_v = v_tail[0], v_tail[1], v_tail[2]

    l_h, c_root_h, c_tip_h, b_h, S_h, x_le_h = _calc_h_tail(xcg_emp, xcg_aft, MAC, S, A_h, l_fuselage, taper_h,
                                                            V_h_norm,
                                                            QC_sweep_h)
    l_v, c_root_v, c_tip_v, b_v, S_v, x_le_v = calc_v_tail(xcg_emp, xcg_aft, b, S, A_v, l_fuselage, taper_v, V_v_norm,
                                                           QC_sweep_v)

    cg_locations = np.array([xcg_fuse, xcg_emp, xcg_fix, xcg_nac, xcg_prop, xcg_wing, xcg_fwd, xcg_aft, zcg])
    tail_h = np.array([l_h, c_root_h, c_tip_h, b_h, S_h])
    tail_v = np.array([l_v, c_root_v, c_tip_v, b_v, S_v])

    X_LE_root = X_LEMAC - Y_MAC * np.tan(LE_sweep)
    alv_h, alv_v = x_le_h - X_LE_root, x_le_v - X_LE_root
    # print("Horizontal tail arm = ", l_h)
    # print("Vertical tail arm = ", l_v)
    # print("Horizontal tail surface fraction = ", (S_h/S))
    # print("Vertical tail surface fraction = ", (S_v/S))
    # print("Fuselage length = ", l_fuselage)
    return cg_locations, tail_h, tail_v, X_LEMAC, alv_h, alv_v


def _calc_h_tail(x_h, xcg_aft, MAC, S, A_h, l_fuselage, tap_h, V_h_norm, sweep_quartchord_h):
    n = 1.
    C_t_h, C_r_h, b_h, S_h = 0, 0, 0, 0
    # print(A_h)
    # print(S)
    # print(V_h_norm)
    while n > 0.001:
        S_frac_h = (V_h_norm * MAC) / (x_h - xcg_aft)
        S_h = S_frac_h * S
        b_h = np.sqrt(A_h * S_h)
        C_r_h = (2. * S_h) / (b_h * (1. + tap_h))
        C_t_h = C_r_h * tap_h
        MAC_h = (2. / 3.) * C_r_h * ((1. + tap_h + tap_h ** 2) / (1. + tap_h))
        sweep_LE_h = np.arctan(np.tan(sweep_quartchord_h) - (C_r_h / (2. * b_h)) * (tap_h - 1))
        sweep_TE_h = np.arctan((C_t_h - C_r_h) / (b_h / 2.) + np.tan(sweep_LE_h))
        y_MAC_h = (b_h / 6.) * ((1. + 2. * tap_h) / (1. + tap_h))
        fuse_MAC_h = MAC_h - 0.25 * MAC_h - y_MAC_h * np.tan(sweep_TE_h)

        n = abs((x_h - (l_fuselage - fuse_MAC_h)) / x_h)
        x_h = l_fuselage - fuse_MAC_h

    x_le_h = l_fuselage - C_r_h
    l_h = x_h - xcg_aft
    
    return l_h, C_r_h, C_t_h, b_h, S_h, x_le_h


def calc_v_tail(x_v, xcg_aft, b, S, A_v, l_fuselage, tap_v, V_v_norm, sweep_quartchord_v):
    n = 1.
    C_r_v, C_t_v, b_v, S_v = 0, 0, 0, 0

    while n > 0.001:
        S_frac_v = (V_v_norm * b) / (x_v - xcg_aft)
        S_v = S_frac_v * S
        b_v = np.sqrt(A_v * S_v)

        C_r_v = S_v / (b_v * (1. + tap_v))
        C_t_v = C_r_v * tap_v

        MAC_v = (2. / 3.) * C_r_v * ((1. + tap_v + tap_v ** 2) / (1. + tap_v))
        sweep_LE_v = np.arctan(np.tan(sweep_quartchord_v) - (C_r_v / (2. * b_v)) * (tap_v - 1))
        y_MAC_v = (b_v / 6.) * ((1. + 2. * tap_v) / (1. + tap_v))

        fuse_MAC_y_new = l_fuselage - C_r_v + 0.25 * MAC_v + np.tan(sweep_LE_v) * y_MAC_v

        sweep_TE_v = np.arctan((C_t_v - C_r_v) / b_v + np.tan(sweep_LE_v))
        fuse_MAC_v = MAC_v - 0.25 * MAC_v - y_MAC_v * np.tan(sweep_TE_v)

        n = abs((x_v - (l_fuselage - fuse_MAC_v)) / (x_v))
        x_v = l_fuselage - fuse_MAC_v

    l_v = x_v - xcg_aft
    x_le_v = l_fuselage - C_r_v

    return l_v, C_r_v, C_t_v, b_v, S_v, x_le_v

# class_I_empennage(MAC, l_fuselage, x_eng, l_n, xcg_OEW_MAC, mass_frac_OEW, xcg_payload, mass_frac_payload, xcg_fuel,
#                  mass_frac_fuel, D_fuse, b)


# def size_tail(wing_area, volume_fraction, tail_sweep, aspect_ratio):
#     tail_area = volume_fraction * wing_area
#     taper_ratio = 0.2 * (2 - tail_sweep * np.pi / 180)
#     tail_span = np.sqrt(tail_area * aspect_ratio)
#     root_chord = (2 * tail_area) / ((1 + taper_ratio) * tail_span)
#     tip_chord = root_chord * taper_ratio
#     return root_chord, tip_chord, tail_span, tail_area
#
#
# def classI_landinggear(MTOW, MLW, l_nosecone):
#     N_mw = MLW / 210000
#     N_mw = int(N_mw)
#
#     while N_mw % 4 != 0:
#         N_mw = N_mw + 1
#
#     if N_mw <= 12:
#         N_s = 2
#     else:
#         N_s = 4
#
#     # l_nlg =
#
#     return N_mw, N_s

# print (classI_landinggear(MTOW, MLW, l_nosecone))
