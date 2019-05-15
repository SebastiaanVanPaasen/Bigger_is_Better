from class_I.class_I_weight_estimation import class_I
from class_I.fuselage_cross_section import fuselage_cross_section
from class_I.planform import wing_parameters, determine_half_chord_sweep
from class_I.drag_estimation import Wing_wetted_area, H_tail_wetted_area, V_tail_wetted_area, Fus_wetted_area, \
    Zero_Lift_Drag_est
from class_I.class_I_empennage_landinggear import class_I_empennage
from class_I.flight_envelope import manoeuvring_envelope, gust_envelope
from avl.conv_wing_avl import make_avl_file, run_avl, find_clalpha
from class_II_weight_estimation import *
from input_files.high_bypass_ratio import *

# M_cruise_list = np.arange(0.5, 0.8, 0.01)
# h_cruise_list = np.arange(3000, 13000, 1000)
# fuel_consumption = np.arange(0.4, 0.9, 0.1)
# aspect_ratios = np.arange(5, 16, 0.5)
# result_wing = []
# result_empty = []
# result_fuel = []
# result_velocity = []
# result_drag = []
# result_cd = []
# result_cl = []

# CL_cruise = np.sqrt((CD_0 * np.pi * A * Oswald) / 3)

# Temp_cruise = Temp_0 + a * h_cruise  # K  based on the altitude you fly at
# a_cruise = np.sqrt(gamma * R_gas * Temp_cruise)  # m/s based on the temperature
#
# Rho_Cruise = Rho_0 * (
#         (1 + (a * h_cruise) / Temp_0) ** (-(g_0 / (R_gas * a))))  # kg/m^3   based on cruise altitude

# V_cruise = M_cruise * a_cruise
# Parameters that will be iterated through-out the process ---------------------------------------------------------
T_input = 0.3
S_input = 6000

# Starting the iteration process -----------------------------------------------------------------------------------
i = 0
maximum = 50
percentage = 10
empty_weight = np.array([])
iteration = {}
total = {}
# final_diagram(CD_0, Oswald)
while i < maximum and percentage > 0.0005:
    print("Starting on iteration: " + str(i))
    # Performing class I weight estimation -------------------------------------------------------------------------
    weights = class_I(CL_cruise, CD_cruise, mission_range, reserve_range, V_cruise, c_j_cruise, W_tfo_frac,
                      W_e_frac,
                      fuel_fractions, N_crew, N_pas, W_person)

    W_TO, W_E_I, W_P, W_F = weights[0], weights[1], weights[2], weights[3]  # N
    iteration["weights"] = [W_TO, W_E_I, W_P, W_F]
    mass_fractions[6] = (weights[1] / weights[0])  # empty mass fraction
    mass_fractions[7] = (weights[2] / weights[0])  # payload mass fraction
    mass_fractions[8] = (weights[3] / weights[0])  # fuel mass fraction

    # Choosing a design point based on the T/W-W/S diagram ---------------------------------------------------------
    T, S = W_TO * T_input, W_TO / S_input
    CL_cruise = W_TO / (0.5 * Rho_Cruise * (V_cruise ** 2) * S)
    # print(Rho_Cruise, V_cruise, S, CD_cruise)
    D = 0.5 * Rho_Cruise * (V_cruise ** 2) * S * CD_cruise
    iteration["thrust, area, cl"] = [T, S, CL_cruise]

    # Sizing the fuselage based on statistical estimated values ----------------------------------------------------
    fuselage_design = fuselage_cross_section(N_pas, N_pas_below)
    iteration["fuselage"] = fuselage_design

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

    # Sizing the wing based on cruise parameters and wing configuration --------------------------------------------
    QC_sweep, LE_sweep, b, c_root, c_tip, dihedral, t_over_c, mac, taper = wing_parameters(M_cruise, CL_cruise, S,
                                                                                           A, wing_option)

    # Calculate additional required half chord sweep for later use
    HC_sweep = determine_half_chord_sweep(c_tip, QC_sweep, c_root, b)
    iteration["planform"] = [QC_sweep, LE_sweep, b, c_root, c_tip, dihedral, t_over_c, mac, HC_sweep]

    # Perform first order cg-range estimation based on statistics --------------------------------------------------
    # Note that the cg-location estimates should be updated after the first iteration!
    x_payload = 0.5 * l_fuselage  # m     cg-location payload w.r.t. nose
    cg_locations, tail_h, tail_v, x_lemac, avl_h, avl_v = class_I_empennage(mass_fractions, mac, l_fuselage,
                                                                            x_engines,
                                                                            l_nacelle, xcg_oew_mac, x_payload,
                                                                            x_fuel,
                                                                            d_fuselage, b, S, taper, v_tail,
                                                                            LE_sweep, h_tail)

    l_h, c_root_h, c_tip_h, b_h, S_h = tail_h[0], tail_h[1], tail_h[2], tail_h[3], tail_h[4]
    l_v, c_root_v, c_tip_v, b_v, S_v = tail_v[0], tail_v[1], tail_v[2], tail_v[3], tail_v[4]

    iteration["cg's"] = [cg_locations, x_lemac]

    # Calculate the accompanying tail sizes ------------------------------------------------------------------------
    HC_sweep_h = determine_half_chord_sweep(c_tip_h, QC_sweep_h, c_root_h, b_h)
    HC_sweep_v = determine_half_chord_sweep(c_tip_v, QC_sweep_v, c_root_v, b_v)

    iteration["horizontal tail"] = [c_root_h, c_tip_h, b_h, S_h, tap_h]
    iteration["vertical tail"] = [c_root_v, c_tip_v, b_v, S_v, tap_v]

    # Calculate new wetted area-------------------------------------------------------------------------------------
    main_wing_wet = Wing_wetted_area(c_root, c_tip, d_fuselage, b, S, winglet_height)
    horizontal_tail_wet = H_tail_wetted_area(c_root_h, tap_h, b_h)
    vertical_tail_wet = V_tail_wetted_area(c_root_v, tap_v, b_v)
    fuselage_wet = Fus_wetted_area(d_fuselage, l_nosecone, l_cabin, l_tailcone)

    # Determine the new CD_0 value ---------------------------------------------------------------------------------
    CD_0 = Zero_Lift_Drag_est(S, main_wing_wet, horizontal_tail_wet, vertical_tail_wet, fuselage_wet)
    # print("the cd_0")
    # print(CD_0)

    # Use AVL to determine the CL_alpha of the wing based on the wing geometry -------------------------------------
    make_avl_file(c_root, c_tip, b, LE_sweep, dihedral, S, CD_0, M_cruise, avl_h, b_h, c_root_h, QC_sweep_h, tap_h,
                  avl_v,
                  b_v, c_root_v, QC_sweep_v, tap_v, tail_type, 12, 5)

    Oswald, CD_cruise = run_avl(CL_cruise, M_cruise, CD_0)
    CL_alpha = (find_clalpha(M_cruise, CD_0) * 180) / np.pi

    iteration["updated parameters"] = [CD_0, Oswald, CL_alpha]

    # Determine maximum loads based on manoeuvring and gust envelope------------------------------------------------
    manoeuvring_loads = manoeuvring_envelope(W_TO, h_cruise, CL_Cruise_max, S, V_cruise)
    gust_loads = gust_envelope(W_TO, h_cruise, CL_alpha, S, mac, V_cruise, manoeuvring_loads[4])

    V_D = manoeuvring_loads[4][3]

    n_max_manoeuvring = max(manoeuvring_loads[1])
    n_max_gust = max(gust_loads[1])

    n_ult = 1.5 * max(n_max_manoeuvring, n_max_gust)

    # Note that for a conventional tail the horizontal tail starts at the root of the vertical tail, thus z_h = 0
    if tail_type == 0:
        z_h = 0
    else:
        z_h = 2

    iteration["n_ult"] = n_ult

    # Start the class II weight estimation -------------------------------------------------------------------------
    # Determine the structural weight components -------------------------------------------------------------------

    t_max_root = t_over_c * c_root
    w_weight = wing_weight(W_TO, W_F, b, HC_sweep, n_ult, S, t_max_root, wing_choice)
    mass_fractions[0] = (w_weight * lbs_to_kg * g_0) / W_TO

    emp_weight = empennage_weight(empennage_choice, np.array([S_h, S_v]), V_D, np.array([HC_sweep_h, HC_sweep_v]),
                                  z_h, b_v)

    mass_fractions[1] = (emp_weight * lbs_to_kg * g_0) / W_TO

    # Note that the fuselage is assumed circular in this case
    fus_weight = fuselage_weight(fuselage_choice, V_D, l_h, d_fuselage, d_fuselage, fuselage_wet)
    mass_fractions[2] = (fus_weight * lbs_to_kg * g_0) / W_TO

    # Note that the choice includes the engine choice here
    nac_weight = nacelle_weight(W_TO, nacelle_choice)
    mass_fractions[3] = (nac_weight * lbs_to_kg * g_0) / W_TO

    lg_weight = landing_gear_weight(W_TO)
    eng_weight = engine_weight(N_engines, w_engine)

    structural_weight = w_weight + emp_weight + fus_weight + nac_weight + lg_weight

    # Determine the propulsion system weight components ------------------------------------------------------------
    ai_weight = induction_weight(duct_length, n_inlets, a_inlets, induction_choice)

    n_prop = prop_characteristics[0]
    n_blades = prop_characteristics[1]
    d_prop = prop_characteristics[2]

    if propeller_choice == 1:
        to_power = ((T * V_cruise) / (550 * prop_characteristics[3]))
        prop_weight = propeller_weight(prop_choice, n_prop, d_prop, to_power, n_blades)
    else:
        prop_weight = 0
        to_power = 0

    # Note that choice is regarding type of fuel tanks
    fuel_sys_weight = fuel_system_weight(N_engines, n_fuel_tanks, W_F, fuel_sys_choice)
    # Choice is depending on type of engine controls and whether there is an afterburner
    w_ec = calc_w_ec(l_fuselage, N_engines, b, engine_choice)
    # Choice is depending on type of starting system and type of engine
    w_ess = calc_w_ess(w_engine, N_engines, start_up_choice)
    # Choice is depending on type of engine
    w_pc = calc_w_pc(n_blades, n_prop, d_prop, N_engines, to_power, prop_choice)
    # Choice is depending on type of engines
    w_osc = calc_w_osc(oil_choice, w_engine, N_engines)

    prop_sys_weight = (w_engine * N_engines) / lbs_to_kg + ai_weight + prop_weight + fuel_sys_weight + w_ec + w_ess \
                      + w_pc + w_osc
    mass_fractions[4] = (prop_sys_weight * lbs_to_kg * g_0) / W_TO

    # Determine fixed equipment weight components ------------------------------------------------------------------
    # Calculate dynamic pressure at dive speed
    q_D = 0.5 * Rho_Cruise * V_D
    w_fc = calc_w_fc(W_TO, q_D)

    # Choice depends on type of engines
    w_hps_els = calc_w_hps_els(hydro_choice, W_TO, V_pax)

    # Maximum range has still to be determined
    w_instr = calc_w_instr(W_E_I, maximum_range)

    w_api = calc_w_api(V_pax, N_crew, N_pas)
    w_ox = calc_w_ox(N_crew, N_pas)
    w_apu = calc_w_apu(W_TO)
    w_fur = calc_w_fur(W_TO, W_F)
    w_bc = calc_w_bc(S_ff)

    fix_equip_weight = w_fc + w_hps_els + w_instr + w_api + w_ox + w_apu + w_fur + w_bc
    mass_fractions[5] = (fix_equip_weight * lbs_to_kg * g_0) / W_TO

    # Determine final operational empty weight ---------------------------------------------------------------------
    W_E_II = (structural_weight + prop_sys_weight + fix_equip_weight) * lbs_to_kg * g_0

    iteration["new empty weight"] = W_E_II
    W_e_frac = W_E_II / W_TO

    percentage = abs((W_E_II - W_E_I) / W_E_I)

    print("the percentage equals: " + str(percentage))
    total[str(i)] = iteration
    i += 1

print("The take-off weight equals: " + str(round(weights[0], 2)))
print("The fuel weight equals: " + str(round(weights[3], 2)))
print("The payload weight equals: " + str(round(weights[2], 2)))
print("The empty weight equals: " + str(round(weights[1], 2)))
print("The new approximated CL_cruise is: " + str(CL_cruise))
print("The new approximated CD_0 is: " + str(CD_0))
print("The new approximated Lift over drag is: " + str(CL_cruise/CD_cruise))
print("The span equals: " + str(b))

# print("The required thrust equals: " + str(T))
# print("The calculated surface area equals: " + str(S))
#
# print("The wing wetted area is: " + str(main_wing_wet))
# print("The vertical tail wetted area is: " + str(vertical_tail_wet))
# print("The horizontal tail wetted area is: " + str(horizontal_tail_wet))
# print("The fuselage wetted area is: " + str(fuselage_wet))
#
# print("The CL_alpha value equals: " + str(CL_alpha))
# print("The new operating empty weight equals: " + str(W_E))

# iterations = np.arange(0, len(empty_weight), 1)
# final_diagram(CD_0, Oswald)
# plt.plot(iterations, empty_weight)
# plt.show()

# result_fuel.append(W_F)
# result_wing.append(w_weight * lbs_to_kg * g_0)
# result_empty.append(W_E_II)
# result_velocity.append(V_cruise)
# result_drag.append(D)
# result_cd.append(CD_cruise)
# result_cl.append(CL_cruise)
#
# plt.figure(1)
# plt.plot(aspect_ratios, result_fuel, label="fuel")
# plt.plot(aspect_ratios, result_wing, label="wing")
# plt.plot(aspect_ratios, result_empty, label="empty")
# plt.plot(h_cruise_list, result_drag, label="drag")
# plt.plot(h_cruise_list, result_velocity, label="velocity")
# plt.plot(h_cruise_list, result_cd, label="cd")
# plt.plot(result_cd, result_cl)
# plt.legend()
# plt.show()
