import numpy as np  # delete later should be loaded in the general program
import matplotlib.pyplot as plt

# option = 0  # 0 for low wing 1 for high wing


def wing_parameters(m_cruise, cl_cruise, surface_area, aspect_ratio, option):
    # Sweep calculation
    m_t = 0.935  # Technology factor for super critical airfoil
    m_dd = m_cruise + 0.03  # drag divergence mach number

    # Torenbeek estimation for quarter chord sweep values are in radians
    if m_cruise <= 0.7:
        quarter_chord_sweep = 0.
    else:
        quarter_chord_sweep = np.arccos(0.75 * (m_t / m_dd))

    # calculate taper ratio using Torenbeek
    taper_ratio = 0.2 * (2. - quarter_chord_sweep)

    # calculate span, root and tip chords geometrically
    span = np.sqrt(surface_area * aspect_ratio)
    chord_root = (2 * surface_area) / ((1 + taper_ratio) * span)
    chord_tip = taper_ratio * chord_root

    # Select dihedral using ADSEE approach
    if option == 0:
        dihedral = np.deg2rad(
            3 - (np.rad2deg(quarter_chord_sweep) / 10) + 2)  # selected for low wing if height wing use -2 instead of +2
    else:
        dihedral = np.deg2rad(
            3 - (np.rad2deg(quarter_chord_sweep) / 10) - 2)

    # thickness over chord ratio
    leading_edge_sweep = np.arctan(np.tan(quarter_chord_sweep) - (chord_root / (2 * span) * (taper_ratio - 1)))

    # print("The sweep at the leading edge equals: " +str(leading_edge_sweep))

    half_chord_sweep = np.arctan(
        ((chord_tip / 2 + span / 2 * np.tan(leading_edge_sweep)) - (chord_root / 2)) / (span / 2))

    thickness_over_chord = (np.cos(half_chord_sweep) ** 3 * (
            m_t - m_dd * np.cos(half_chord_sweep)) - 0.115 * cl_cruise ** 1.5) / (np.cos(half_chord_sweep) ** 2)

    if thickness_over_chord < 0.10:
        thickness_over_chord = 0.10
    else:
        thickness_over_chord = min(thickness_over_chord, 0.18)

    mac = (2 / 3) * chord_root * ((1 + taper_ratio + (taper_ratio ** 2)) / (1 + taper_ratio))

    # print("here again")
    # print(leading_edge_sweep)
    # print(m_t - m_dd)
    # print(half_chord_sweep)
    # print(cl_cruise)
    # print(thickness_over_chord)

    return (
        quarter_chord_sweep, leading_edge_sweep, span, chord_root, chord_tip, dihedral,
        thickness_over_chord, mac, taper_ratio)
    
# def wing_parameters_boxed(m_cruise, cl_cruise, surface_area, aspect_ratio):
        
    




# inputs
# s_ratio = 0.05
# v_tail_le_sweep = np.rad2deg(35)
# v_tail_aspect_ratio = 1.5
# v_tail_taper_ratio = 0.4
# h_tail_height = 0.5 #ratio of tail height on vertical tail 0.5 is at half the vertical tail
#    

def tail_distance(s_ref, s_ratio, aspect_ratio, leadin_edge_sweep, h_tail_height):
    s = s_ref * s_ratio
    span = np.sqrt(s * aspect_ratio)
    height = h_tail_height * span
    x_dist_tails = height / np.tan(leadin_edge_sweep)
    # chord_root = (2 * s) / ((1 + taper_ratio) * span)
    # chord_tip = taper_ratio * chord_root
    # print(chord_root, chord_tip, x_dist_tails)
    return (x_dist_tails)


def plot_planform(leading_edge_sweep, chord_root, chord_tip, span):
    x_list = [0, span / 2, span / 2, 0, 0]
    y_list = [0, span / 2 * np.tan(leading_edge_sweep), span / 2 * np.tan(leading_edge_sweep) + chord_tip, chord_root,
              0]
    plt.plot(x_list, y_list)
    plt.axis([0, 35, 0, 30], 'equal')

    plt.plot([0, chord_root], [0, 0])
    plt.plot([span / 2 * np.tan(leading_edge_sweep), span / 2 * np.tan(leading_edge_sweep) + chord_tip],
             [span / 2, span / 2])
    plt.plot([0, span / 2 * np.tan(leading_edge_sweep)], [0, span / 2])
    plt.plot([chord_root, span / 2 * np.tan(leading_edge_sweep) + chord_tip], [0, span / 2])
    plt.show()


def determine_half_chord_sweep(chord_tip, qc_sweep, chord_root, span):
    hc_sweep = np.arctan(((chord_tip / 4 + span / 2 * np.tan(qc_sweep)) - (chord_root / 4)) / (span / 2))

    return hc_sweep

# M_cruise = 0.7  # inputs from different part
# surface_area = 350  # inputs from different part
# aspect_ratio = 15  # inputs from different part
# CL_cruise = 0.7
#
# quarter_cord_sweep, leading_edge_sweep, taper_ratio, span, cord_root, cord_tip, dihedral, tickness_over_cord, mac = wing_parameters(
#     M_cruise, CL_cruise, surface_area, aspect_ratio)
##print(quarter_cord_sweep, leading_edge_sweep, taper_ratio, span, cord_root, cord_tip, dihedral, tickness_over_cord, mac)
# x_dist_tails = tail_distance(surface_area, s_ratio, v_tail_aspect_ratio, v_tail_le_sweep,  h_tail_height)
# print(x_dist_tails)
