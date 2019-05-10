from inputs import *
from class_I_weight_estimation import class_I
from fuselage_cross_section import fuselage_cross_section
from planform import wing_parameters, determine_half_chord_sweep
from drag_estimation import Wing_wetted_area, H_tail_wetted_area, V_tail_wetted_area, Fus_wetted_area, \
    Zero_Lift_Drag_est
from class_I_empennage_landinggear import class_I_empennage, size_tail
from flight_envelope import manoeuvring_envelope, gust_envelope
from class_II_weight_estimation import *

# parameters that will be updated CD_0, Oswald, W_e_frac, All weights, t_over_c,


# Performing class I weight estimation ---------------------------------------------------------------------------------
weights = class_I(A, Oswald, CD_0, mission_range, reserve_range, V_cruise, c_j_cruise, W_tfo_frac, W_e_frac,
                  fractions, N_crew, N_pas, W_person)

W_TO, W_E, W_P, W_F = weights[0], weights[1], weights[2], weights[3]
m_e_frac = (weights[1] / weights[0])  # empty mass fraction
m_p_frac = (weights[2] / weights[0])  # payload mass fraction
m_f_frac = (weights[3] / weights[0])  # fuel mass fraction

# print("The take-off mass equals: " + str(round(masses[0], 2)))
# print("The fuel mass equals: " + str(round(masses[3], 2)))
# print("The payload mass equals: " + str(round(masses[2], 2)))
# print("The empty mass equals: " + str(round(masses[1], 2)))

# Choosing a design point based on the T/W-W/S diagram -----------------------------------------------------------------
T, S = W_TO * T_input, W_TO / S_input
CL_cruise = W_TO / (0.5 * Rho_Cruise * (V_cruise ** 2) * S)
# print("The required thrust equals: " + str(T))
# print("The calculated surface area equals: " + str(S))
# print("The new approximated CL_cruise is: " + str(CL_cruise))

# Sizing the fuselage based on statistical estimated values ------------------------------------------------------------
fuselage_design = fuselage_cross_section(N_pas, N_pas_below)

# Define fuselage parameters out of the fuselage design
d_fuselage = fuselage_design[1]
l_nosecone = np.sum(fuselage_design[6]) / 2
l_cabin = fuselage_design[2]
l_tailcone = np.sum(fuselage_design[5]) / 2
l_fuselage = fuselage_design[7]

# Define fineness ratio's based on the dimensions that came out of the fuselage design
f = l_fuselage / d_fuselage
f_tc = l_tailcone / d_fuselage
f_nc = l_nosecone / d_fuselage
V_pax = 0.25 * np.pi * (d_fuselage ** 2) * l_cabin  # m^3  Volume of the passenger cabin

# Area of fuselage without cut-outs and wing
# S_fus_gross = np.pi * d_fuselage * l_fuselage * (1 - f_nc / (3 * f) - f_tc / (2 * f))
# Area of freight floor, very rough initial guess for now
S_ff = 0.5 * d_fuselage * l_cabin

# Sizing the wing based on cruise parameters and wing configuration ----------------------------------------------------
QC_sweep, b, c_root, c_tip, dihedral, t_over_c, mac = wing_parameters(M_cruise, CL_cruise, S, A, wing_option)

# Calculate additional required half chord sweep for later use
HC_sweep = determine_half_chord_sweep(c_root, c_tip, QC_sweep, b)

# Perform first order cg-range estimation based on statistics ----------------------------------------------------------
# Note that the cg-location estimates should be updated after the first iteration!
x_payload = 0.5 * l_fuselage  # m     cg-location payload w.r.t. nose
cg_locations, S_h_frac, S_v_frac, x_lemac, l_h = class_I_empennage(mac, l_fuselage, x_engines, l_nacelle, xcg_oew_mac,
                                                                   m_e_frac, x_payload, m_p_frac, x_fuel, m_f_frac,
                                                                   d_fuselage, b)
# print("The resulting geometric locations are: " + str(cg_locations))

# Calculate the accompanying tail sizes --------------------------------------------------------------------------------
c_root_h, c_tip_h, b_h, S_h = size_tail(S, S_h_frac, QC_sweep_h, A_h)
c_root_v, c_tip_v, b_v, S_v = size_tail(S, S_v_frac, QC_sweep_v, A_v)
HC_sweep_h = determine_half_chord_sweep(c_root_h, c_tip_h, b_h, QC_sweep_h)
HC_sweep_v = determine_half_chord_sweep(c_root_v, c_tip_v, b_v, QC_sweep_v)

taper_h = c_root_h / c_tip_h
taper_v = c_root_v / c_tip_v

# Calculate new wetted area---------------------------------------------------------------------------------------------
main_wing_wet = Wing_wetted_area(c_root, c_tip, d_fuselage, b, S, winglet_height)
horizontal_tail_wet = H_tail_wetted_area(c_root_h, taper_h, b_h)
vertical_tail_wet = V_tail_wetted_area(c_root_v, taper_v, b_v)
fuselage_wet = Fus_wetted_area(d_fuselage, l_nosecone, l_cabin, l_tailcone)

# print("The wing wetted area is: " + str(main_wing_wet))
# print("The vertical tail wetted area is: " + str(vertical_tail_wet))
# print("The horizontal tail wetted area is: " + str(horizontal_tail_wet))
# print("The fuselage wetted area is: " + str(fuselage_wet))

# Determine the new CD_0 value -----------------------------------------------------------------------------------------
CD_0 = Zero_Lift_Drag_est(S, main_wing_wet, horizontal_tail_wet, vertical_tail_wet, fuselage_wet)
# print("The new approximated CD_0 is: " + str(CD_0))

# lift_capability = Clean_Wing_Lift(A, M_cruise, HC_sweep, QC_sweep, alpha, alpha_0, CL_des, airfoil_efficiency)
# CL_alpha = lift_capability[0]

# Determine maximum loads based on manoeuvring and gust envelope--------------------------------------------------------
manoeuvring_loads = manoeuvring_envelope(W_TO, h_cruise, CL_Cruise_max, S, V_cruise)
gust_loads = gust_envelope(W_TO, h_cruise, CL_alpha, S, mac, V_cruise, manoeuvring_loads[4])

V_D = manoeuvring_loads[4][3]

n_max_manoeuvring = max(manoeuvring_loads[1])
n_max_gust = max(gust_loads[1])

n_ult = 1.5 * max(n_max_manoeuvring, n_max_gust)

# Note that this is for a conventional tail
z_h = 0

# Start the class II weight estimation ---------------------------------------------------------------------------------
# Determine the structural weight components ---------------------------------------------------------------------------
t_max_root = t_over_c * c_root

w_weight = wing_weight(W_TO, W_F, b, HC_sweep, n_ult, S, t_max_root, [0, 0, 0, 0, 0, 0])
emp_weight = empennage_weight(np.array([0, 0]), np.array([S_h, S_v]), V_D, np.array([HC_sweep_h, HC_sweep_v]), z_h, b_v)

# Note that the fuselage is assumed circular in this case
fus_weight = fuselage_weight(0, V_D, l_h, d_fuselage, d_fuselage, fuselage_wet)

# Note that the choice includes the engine choice here
nac_weight = nacelle_weight(W_TO, 0)
lg_weight = landing_gear_weight(W_TO)
eng_weight = engine_weight(N_engines, w_engine)

structural_weight = w_weight + emp_weight + fus_weight + nac_weight + lg_weight + eng_weight

# Determine the propulsion system weight components --------------------------------------------------------------------
ai_weight = induction_weight(duct_length, n_inlets, a_inlets, [0, 0])

n_prop = prop_characteristics[0]
n_blades = prop_characteristics[1]
d_prop = prop_characteristics[2]

if propeller_choice == 1:
    to_power = ((T * V_cruise) / (550 * prop_characteristics[3])) / hp_to_N
    prop_weight = propeller_weight(0, n_prop, d_prop, to_power, n_blades)
else:
    prop_weight = 0
    to_power = 0

# Note that choice is regarding type of fuel tanks
fuel_sys_weight = fuel_system_weight(N_engines, n_fuel_tanks, W_F, 0)
# Choice is depending on type of engine controls and whether there is an afterburner
w_ec = calc_w_ec(l_fuselage, N_engines, b, [2, 0])
# Choice is depending on type of starting system and type of engine
w_ess = calc_w_ess(w_engine, N_engines, 1)
# Choice is depending on type of engine
w_pc = calc_w_pc(n_blades, n_prop, d_prop, N_engines, to_power, 1)
# Choice is depending on type of engines
w_osc = calc_w_osc(1, w_engine, N_engines)

prop_sys_weight = w_engine * N_engines + ai_weight + prop_weight + fuel_sys_weight + w_ec + w_ess + w_pc + w_osc

# Determine fixed equipment weight components --------------------------------------------------------------------------
# Calculate dynamic pressure at dive speed
q_D = 0.5 * Rho_Cruise * V_D
w_fc = calc_w_fc(W_TO, q_D)
# Choice depends on type of engines
w_hps_els = calc_w_hps_els(0, W_TO, V_pax)
# Maximum range has still to be determined
w_instr = calc_w_instr(W_E, maximum_range)

w_api = calc_w_api(V_pax, N_crew, N_pas)
w_ox = calc_w_ox(N_crew, N_pas)
w_apu = calc_w_apu(W_TO)
w_fur = calc_w_fur(W_TO, W_F)
w_bc = calc_w_bc(S_ff)

fix_equip_weight = w_fc + w_hps_els + w_instr + w_api + w_ox + w_apu + w_fur + w_bc

# Determine final operational empty weight -----------------------------------------------------------------------------
W_E = structural_weight + prop_sys_weight + fix_equip_weight

print("The new operating empty weight equals: " + str(W_E * lbs_to_kg * g_0))
