from inputs import *
from class_I_weight_estimation import class_I
from fuselage_cross_section import fuselage_cross_section
from planform import wing_parameters
from drag_estimation import Wing_wetted_area, H_tail_wetted_area, V_tail_wetted_area, Fus_wetted_area, \
    Zero_Lift_Drag_est
from lift_estimation import Clean_Wing_Lift

weights = class_I(A, Oswald, S_ratio, C_fe, cruise_range, reserve_range, V_cruise, c_j_cruise, W_tfo_frac, W_e_frac,
                  fractions, N_crew, N_pas, W_person)
# print("The take-off mass equals " + str(round(weights[0], 2)))
# print("The fuel mass equals " + str(round(weights[3], 2)))
# print("The payload mass equals " + str(round(weights[2], 2)))
# print("The empty mass equals " + str(round(weights[1], 2)))

design_point = np.array([0.3, 6000])

T, S = weights[0] * design_point[0], weights[0] / design_point[1]

fuselage_design = fuselage_cross_section(N_pas, N_pas)
QC_sweep, LE_sweep, taper_ratio, b, root_chord, tip_chord, dihedral, t_over_c = wing_parameters(M_cruise, CL_cruise, S, A)

d_fuselage = fuselage_design[1]
l_nosecone = np.sum(fuselage_design[5]) / 2
l_cabin = fuselage_design[2]
l_tailcone = np.sum(fuselage_design[6]) / 2

winglet_height = 0.
# change function and add winglet height for the boxed wing
wing_wetted_area = Wing_wetted_area(root_chord, tip_chord, d_fuselage, b, S, winglet_height)
H_wetted_area = H_tail_wetted_area(root_chord_h, taper_ratio_h, span_h)
V_wetted_area = V_tail_wetted_area(root_chord_v, taper_ratio_v, span_v)
Fus_wetted_area = Fus_wetted_area(d_fuselage, l_nosecone, l_cabin, l_tailcone)
CD_0 = Zero_Lift_Drag_est(S, Wing_wetted_area, H_wetted_area, V_wetted_area, Fus_wetted_area)
lift_capability = Clean_Wing_Lift(A, M_cruise, HC_sweep, QC_sweep, alpha, alpha_0, CL_des, airfoil_efficiency)

# print("The wing wetted area is: " + str(Wing_wetted_area))
# print("The vertical tail wetted area is: " + str(V_wetted_area))
# print("The horizontal tail wetted area is: " + str(H_wetted_area))
# print("The fuselage wetted area is: " + str(Fus_wetted_area))
# print("The new approximated CD_0 is: " + str(CD_0))
# aspect_ratio, m_cruise, half_chord_sweep, quarter_chord_sweep, alpha, alpha_0, cl_des, airfoil_eff

