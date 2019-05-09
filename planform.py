import numpy as np  # delete later should be loaded in the general program
import matplotlib.pyplot as plt


def wing_parameters(m_cruise, cl_cruise, surface_area, aspect_ratio):
    # Sweep calculation
    m_t = 0.935  # Technology factor for super critical airfoil
    m_dd = m_cruise + 0.03  # drag divergence mach number

    # Torenbeek estimation for quarter chord sweep values are in radians
    if m_cruise <= 0.69:
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
    dihedral = np.deg2rad(
        3 - (np.rad2deg(quarter_chord_sweep) / 10) + 2)  # selected for low wing if height wing use -2 instead of +2

    # thickness over chord ratio
    leading_edge_sweep = np.arctan(np.tan(quarter_chord_sweep) - (chord_root / (2 * span) * (taper_ratio - 1)))

    # print("The sweep at the leading edge equals: " +str(leading_edge_sweep))

    half_chord_sweep = np.arctan(
        ((chord_tip / 2 + span / 2 * np.tan(leading_edge_sweep)) - (chord_root / 2)) / (span / 2))
    thickness_over_chord = (np.cos(half_chord_sweep) ** 3 * (
            m_t - m_dd * np.cos(half_chord_sweep)) - 0.115 * cl_cruise ** 1.5) / (np.cos(half_chord_sweep) ** 2)
    thickness_over_chord = min(thickness_over_chord, 0.18)

    return (
        quarter_chord_sweep, leading_edge_sweep, taper_ratio, span, chord_root, chord_tip, dihedral,
        thickness_over_chord)


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


def determine_half_chord_sweep(root_chord, tip_chord, span, qc_sweep):
    root_half = root_chord * 0.5
    tip_half = root_chord * 0.25 * (span / 2) * np.tan(np.radians(qc_sweep)) + 0.25 * tip_chord
    hc_sweep = np.tan((root_half - tip_half) / (span / 2))

# M_cruise = 0.7  # inputs from different part
# surface_area = 350  # inputs from different part
# aspect_ratio = 15  # inputs from different part
# CL_cruise = 0.7
#
# quarter_cord_sweep, leading_edge_sweep, taper_ratio, span, cord_root, cord_tip, dihedral, tickness_over_cord = wing_parameters(
#     M_cruise, CL_cruise, surface_area, aspect_ratio)
# print(quarter_cord_sweep, leading_edge_sweep, taper_ratio, span, cord_root, cord_tip, dihedral, tickness_over_cord)
