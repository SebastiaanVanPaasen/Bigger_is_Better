from inputs import *
from class_I_weight_estimation import class_I
from fuselage_cross_section import fuselage_cross_section
from planform import wing_parameters
from drag_estimation import Wing_wetted_area

weights = class_I(A, Oswald, S_ratio, C_fe, cruise_range, reserve_range, V_cruise, c_j_cruise, W_tfo_frac, W_e_frac,
                  fractions, N_crew, N_pas, W_person)
# print("The take-off mass equals " + str(round(weights[0], 2)))
# print("The fuel mass equals " + str(round(weights[3], 2)))
# print("The payload mass equals " + str(round(weights[2], 2)))
# print("The empty mass equals " + str(round(weights[1], 2)))

design_point = np.array([0.3, 6000])

T, S = weights[0] * design_point[0], weights[0] / design_point[1]

fuselage_design = fuselage_cross_section(N_pas, N_pas)
QC_sweep, LE_sweep, taper_ratio, b, c_root, c_tip, dihedral, t_over_c = wing_parameters(M_cruise, CL_cruise, S, A)

d_fuselage = fuselage_design[1]
l_nosecone = np.sum(fuselage_design[5]) / 2
l_tailcone = np.sum(fuselage_design[6]) / 2

winglet_height = 0.
wing_wetted_area = Wing_wetted_area(c_root, c_tip, d_fuselage, b, S,
                                    winglet_height)  # change function and add winglet height for the boxed wing

L1 = 5
L2 = 30
L3 = 5

half_chord_sweep = 0.5
quarter_chord_sweep = 0.6
alpha = 0.5
alpha_0 = 0.3
Cl_des = 1.2
airfoil_eff_factor = 0.95

span = np.sqrt(S * A)
root_chord_h = 5
root_chord_v = 4
span_h = 20
span_v = 25
taper_ratio_h = 0.5
taper_ratio_v = 0.5
S_h = 50
S_v = 70

H_wetted_area = H_tail_wetted_area(root_chord_h, taper_ratio_h, span_h)
V_wetted_area = V_tail_wetted_area(root_chord_v, taper_ratio_v, span_v)
Fus_wetted_area = Fus_wetted_area(D_fuselage, L1, L2, L3)
CD_0 = Zero_Lift_Drag_est(S, Wing_wetted_area, H_wetted_area, V_wetted_area, Fus_wetted_area)

print(Wing_wetted_area, H_wetted_area, V_wetted_area, Fus_wetted_area, CD_0)
