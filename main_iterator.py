from inputs import *
from class_I_weight_estimation import class_I
from fuselage_cross_section import fuselage_cross_section
from planform import wing_parameters, determine_half_chord_sweep
from drag_estimation import Wing_wetted_area, H_tail_wetted_area, V_tail_wetted_area, Fus_wetted_area, \
    Zero_Lift_Drag_est
from lift_estimation import Clean_Wing_Lift
from class_I_empennage_landinggear import class_I_empennage, size_tail

masses = class_I(A, Oswald, S_ratio, C_fe, cruise_range, reserve_range, V_cruise, c_j_cruise, W_tfo_frac, W_e_frac,
                 fractions, N_crew, N_pas, W_person)
m_e_frac = (masses[1] / masses[0])
m_p_frac = (masses[2] / masses[0])
m_f_frac = (masses[3] / masses[0])

print("The take-off mass equals: " + str(round(masses[0], 2)))
print("The fuel mass equals: " + str(round(masses[3], 2)))
print("The payload mass equals: " + str(round(masses[2], 2)))
print("The empty mass equals: " + str(round(masses[1], 2)))

design_point = np.array([0.3, 6000])
T, S = masses[0] * design_point[0], masses[0] / design_point[1]
print("The required thrust equals: " + str(T))
print("The calculated surface area equals: " + str(S))

fuselage_design = fuselage_cross_section(N_pas, N_pas)
QC_sweep, LE_sweep, taper_ratio, b, root_chord, tip_chord, dihedral, t_over_c = wing_parameters(M_cruise, CL_cruise, S,
                                                                                                A)
HC_sweep = determine_half_chord_sweep(root_chord, tip_chord, QC_sweep, b)

d_fuselage = fuselage_design[1]
l_nosecone = np.sum(fuselage_design[6]) / 2
l_cabin = fuselage_design[2]
l_tailcone = np.sum(fuselage_design[5]) / 2
l_fuselage = fuselage_design[7]

x_payload = 0.5 * l_fuselage  # m     cg-location payload w.r.t. nose
cg_locations, S_h_frac, S_v_frac, x_lemac = class_I_empennage(mac, l_fuselage, x_engines, l_nacelle, xcg_oew_mac,
                                                              m_e_frac, x_payload, m_p_frac, x_fuel, m_f_frac,
                                                              d_fuselage, b)

root_chord_h, tip_chord_h, taper_ratio_h, span_h = size_tail(S, S_h_frac, QC_sweep_h, A_h)
root_chord_v, tip_chord_v, taper_ratio_v, span_v = size_tail(S, S_v_frac, QC_sweep_v, A_v)
print("The resulting geometric locations are: " + str(cg_locations))

winglet_height = 0.
# change function and add winglet height for the boxed wing
main_wing_wet = Wing_wetted_area(root_chord, tip_chord, d_fuselage, b, S, winglet_height)
horizontal_tail_wet = H_tail_wetted_area(root_chord_h, taper_ratio_h, span_h)
vertical_tail_wet = V_tail_wetted_area(root_chord_v, taper_ratio_v, span_v)
fuselage_wet = Fus_wetted_area(d_fuselage, l_nosecone, l_cabin, l_tailcone)
CD_0 = Zero_Lift_Drag_est(S, main_wing_wet, horizontal_tail_wet, vertical_tail_wet, fuselage_wet)
lift_capability = Clean_Wing_Lift(A, M_cruise, HC_sweep, QC_sweep, alpha, alpha_0, CL_des, airfoil_efficiency)

print("The wing wetted area is: " + str(main_wing_wet))
print("The vertical tail wetted area is: " + str(vertical_tail_wet))
print("The horizontal tail wetted area is: " + str(horizontal_tail_wet))
print("The fuselage wetted area is: " + str(fuselage_wet))
print("The new approximated CD_0 is: " + str(CD_0))
