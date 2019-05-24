import numpy as np

# def Exposed_Areas(root_chord, span, S, taper_ratio, wing type): #exposed wing area calculation
#    projected_length_fuselage = D_fuselage/2*(np.tan(LE_sweep)+root_chord*taper_ratio*2/span)
#    total_area = projected_length_fuselage*D_fuselage/2
#    sweep_excluding_areas = .5*(D_fuselage/2)**2*(np.tan(LE_sweep)+projected_length_fuselage/D_fuselage*2-root_chord)
#    exposed_area = S - 2*(total_area - sweep_excluding_areas)
#    if wing type == "Vertical":
#        exposed_area= exposed_area/2
#    return exposed_area

# D_fuselage = 3
# L1 = 5
# L2 = 30
# L3 = 5
# S = 500
# A = 10
# M_cruise = 0.7
# half_chord_sweep = 0.5
# quarter_chord_sweep = 0.6
# LE_sweep = 0.65
# alpha = 0.5
# alpha_0 = 0.3
# Cl_des = 1.2
# airfoil_eff_factor = 0.95
# root_chord = 15
# span = np.sqrt(S * A)
# taper_ratio = 0.3
# root_chord_h = 5
# root_chord_v = 4
# span_h = 20
# span_v = 25
# taper_ratio_h = 0.5
# taper_ratio_v = 0.5
# S_h = 50
# S_v = 70
# tip_chord = root_chord * taper_ratio
# winglet_height = 0.


def Wing_wetted_area(root_chord, tip_chord, D_fuselage, span, S,):  
    # set winglet area to zero if not needed, select correct input for boxed wing
    fuselage_chord = root_chord * (1 - (D_fuselage / span))  # determine the cord at the outer fuselage

    exposed_area = S - (0.5 * (
            fuselage_chord + root_chord)) * D_fuselage  # calulate the exposed area, substracting the area inside the plane from the S
    wetted_area = 1.07 * 2 * exposed_area  # wetted area approximation
    return wetted_area

def Wing_wetted_area_boxed(root_chord, tip_chord, D_fuselage, span, S,
                     winglet_height):  # set winglet area to zero if not needed, select correct input for boxed wing
    root_chord_single = root_chord/2    
    tip_chord_single = tip_chord/2              
    fuselage_chord = root_chord_single * (1 - (D_fuselage / span))  # determine the cord at the outer fuselage
    winglet_area = tip_chord_single * winglet_height
    exposed_area = S - (0.5 * (
            fuselage_chord + root_chord)) * D_fuselage + winglet_area  # calulate the exposed area, substracting the area inside the plane from the S
    wetted_area = 1.07 * 2 * exposed_area  # wetted area approximation
    return wetted_area


# function to calulate horizontal tail wetted area using simple geometry same goes for the vertical tail
def H_tail_wetted_area(root_chord_h, taper_ratio_h, span_h):
    exposed_area = 0.5 * (root_chord_h * (1 + taper_ratio_h)) * span_h
    wetted_area = 1.05 * 2 * exposed_area
    return wetted_area


def V_tail_wetted_area(root_chord_v, taper_ratio_v, span_v):
    exposed_area = 0.5 * (root_chord_v * (1 + taper_ratio_v)) * (span_v / 2)
    wetted_area = 1.05 * 2 * exposed_area
    return wetted_area


# using an adsee relation the fuselage wetted area is calculated
def Fus_wetted_area(D_fuselage, L1, L2, L3):
    wetted_area = np.pi * D_fuselage / 4 * (1 / (3 * L1 ** 2) * ((4 * L1 ** 2 + D_fuselage ** 2 / 4) ** 1.5) \
                                            - D_fuselage + 4 * L2 + 2 * np.sqrt(L3 ** 2 + D_fuselage ** 2 / 4))
    return wetted_area


# zero lift drag is estimated using wetted areas and factors from adsee, lift induced drag is done in AVL
def Zero_Lift_Drag_est(S, wing_wetted_area, h_wetted_area, v_wetted_area, fus_wetted_area):
    CD_0 = 1.15 * (1 / S) * (0.003 * wing_wetted_area + 0.0024 * fus_wetted_area + 0.0025 * (
            h_wetted_area + v_wetted_area))  # adsee relations add 15% instead of 10 because nacelles are neglected
    return CD_0
