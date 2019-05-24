#import matplotlib.pyplot as plt
import input_files.strutted_wing as inputs
import numpy as np
import csv
#import sys
#sys.path.append("C:/Users/sebas/OneDrive/Documents/DSE/Bigger_is_Better/class_I")

from constants_and_conversions import Temp_0, g_0, gamma, a, R_gas, Rho_0, lbs_to_kg
from class_II_weight_estimation import Class_II
from class_I.class_I_weight_estimation import class_I
from class_I.fuselage_cross_section import fuselage_cross_section
from class_I.planform import wing_parameters, determine_half_chord_sweep
from class_I.drag_estimation import Wing_wetted_area, H_tail_wetted_area, V_tail_wetted_area, Fus_wetted_area, \
    Zero_Lift_Drag_est
from class_I.class_I_empennage_landinggear import class_I_empennage
from class_I.flight_envelope import manoeuvring_envelope, gust_envelope
from avl.conv_wing_avl import make_avl_file, run_avl, find_clalpha


#        if h_cruise < 11000.:
#            Temp_cruise = Temp_0 + a * h_cruise  # K  based on the altitude you fly at
#        else:
#            Temp_cruise = 216.65
#
#        a_cruise = np.sqrt(gamma * R_gas * Temp_cruise)  # m/s based on the temperature
#
#        Rho_Cruise = Rho_0 * (
#                (1 + (a * h_cruise) / Temp_0) ** (-(g_0 / (R_gas * a))))  # kg/m^3   based on cruise altitude
#
#        V_cruise = M_cruise * a_cruise
#        Ct0_value = 12e-06

#        print("vcruise" + str(V_cruise))
#        print("rho cruise" + str(Rho_Cruise))
#        print("a cruise" + str(a_cruise))



#writeFile.close()

def main_iterator(CL_cruise, CD_cruise, W_e_frac, fuel_frac, mass_frac, C_t0, V_cruise):
    with open('people.csv', 'w') as writeFile:
        writer = csv.writer(writeFile)
    
        mass_fractions = inputs.mass_fractions_input
        percentages = []
        empty_weight = []
    
    
        fuselage_design = fuselage_cross_section(inputs.N_pas, inputs.N_pas_below)
        
        # Define fuselage parameters out of the fuselage design
        d_fuselage = fuselage_design[1]
        l_nosecone = np.sum(fuselage_design[6]) / 2
        l_cabin = fuselage_design[2]
        l_tailcone = np.sum(fuselage_design[5]) / 2
        l_fuselage = fuselage_design[7]
        V_pax = 0.25 * np.pi * (d_fuselage ** 2) * l_cabin  
        # m^3  Volume of the passenger cabin
        writer.writerow(["fuselage diameter", d_fuselage])
        writeFile.close()
        # Define fineness ratio's based on the fuselage design
    #    f = l_fuselage / d_fuselage
    #    f_tc = l_tailcone / d_fuselage
    #    f_nc = l_nosecone / d_fuselage
    
        # Area of freight floor, very rough initial guess for now
        S_ff = 0.5 * d_fuselage * l_cabin
        
        # Starting the iteration process ------------------------------------------
        i = 0
        maximum = 500
        percentage = 1
        iteration = {}
    
        while i < maximum and percentage > 0.005:
            print("Starting on iteration: " + str(i))
            # Performing class I weight estimation -------------------------------------------------------------------------
            weights = class_I(CL_cruise, CD_cruise, inputs.mission_range, inputs.reserve_range, V_cruise, inputs.c_j_cruise, inputs.W_tfo_frac,
                              W_e_frac,
                              fuel_frac, inputs.N_pas, inputs.N_crew, inputs.W_person, inputs.W_carg)
    
            W_TO, W_E_I, W_P, W_F = weights[0], weights[1], weights[2], weights[3]  # N
            
            # print(W_TO)
            # print(CL_cruise)
            # print(CD_cruise)
            # print(mission_range)
            # print(V_cruise)
            # print(W_tfo_frac)
            
            iteration["weights"] = [W_TO, W_E_I, W_P, W_F]
            mass_fractions[6] = (weights[1] / weights[0])  # empty mass fraction
            mass_fractions[7] = (weights[2] / weights[0])  # payload mass fraction
            mass_fractions[8] = (weights[3] / weights[0])  # fuel mass fraction
    
            # Choosing a design point based on the T/W-W/S diagram ---------------------------------------------------------
            T, S = W_TO * inputs.T_input, W_TO / inputs.S_input
            CL_cruise = (0.75 * W_TO) / (0.5 * inputs.Rho_Cruise_input * (V_cruise ** 2) * S)
            if CL_cruise > 1.:
                CL_cruise = 1.
    
            # print(Rho_Cruise, V_cruise, S, CD_cruise)
    #            D = 0.5 * Rho_Cruise * (V_cruise ** 2) * S * CD_cruise
            iteration["thrust, area, cl"] = [T, S, CL_cruise]
    
            # Sizing the fuselage based on statistical estimated values ----------------------------------------------------
    
            iteration["fuselage"] = fuselage_design
    
    
            # Sizing the wing based on cruise parameters and wing configuration --------------------------------------------
            QC_sweep, LE_sweep, b, c_root, c_tip, dihedral, t_over_c, mac, taper = wing_parameters(inputs.M_cruise, CL_cruise,
                                                                                                   S,
                                                                                                   inputs.A, inputs.wing_option)
    
            # Calculate additional required half chord sweep for later use
            HC_sweep = determine_half_chord_sweep(c_tip, QC_sweep, c_root, b)
            iteration["planform"] = [QC_sweep, LE_sweep, b, c_root, c_tip, dihedral, t_over_c, mac, HC_sweep]
    
            # Perform first order cg-range estimation based on statistics --------------------------------------------------
            # Note that the cg-location estimates should be updated after the first iteration!
            x_payload = 0.5 * l_fuselage  # m     cg-location payload w.r.t. nose
            cg_locations, tail_h, tail_v, x_lemac, avl_h, avl_v = class_I_empennage(mass_fractions, mac, l_fuselage,
                                                                                    inputs.x_engines,
                                                                                    inputs.l_nacelle, inputs.xcg_oew_mac, x_payload,
                                                                                    inputs.x_fuel,
                                                                                    d_fuselage, b, S, taper, inputs.v_tail,
                                                                                    LE_sweep, inputs.h_tail)
    
            l_h, c_root_h, c_tip_h, b_h, S_h = tail_h[0], tail_h[1], tail_h[2], tail_h[3], tail_h[4]
            l_v, c_root_v, c_tip_v, b_v, S_v = tail_v[0], tail_v[1], tail_v[2], tail_v[3], tail_v[4]
            iteration["cg's"] = [cg_locations, x_lemac]
    
            # Calculate the accompanying tail sizes ------------------------------------------------------------------------
            HC_sweep_h = determine_half_chord_sweep(c_tip_h, inputs.QC_sweep_h, c_root_h, b_h)
            HC_sweep_v = determine_half_chord_sweep(c_tip_v, inputs.QC_sweep_v, c_root_v, b_v)
    
            iteration["horizontal tail"] = [c_root_h, c_tip_h, b_h, S_h, inputs.tap_h]
            iteration["vertical tail"] = [c_root_v, c_tip_v, b_v, S_v, inputs.tap_v]
    
            # Calculate new wetted area-------------------------------------------------------------------------------------
            main_wing_wet = Wing_wetted_area(c_root, c_tip, d_fuselage, b, S, inputs.winglet_height)
            horizontal_tail_wet = H_tail_wetted_area(c_root_h, inputs.tap_h, b_h)
            vertical_tail_wet = V_tail_wetted_area(c_root_v, inputs.tap_v, b_v)
            fuselage_wet = Fus_wetted_area(d_fuselage, l_nosecone, l_cabin, l_tailcone)
    
            # Determine the new CD_0 value ---------------------------------------------------------------------------------
            CD_0 = Zero_Lift_Drag_est(S, main_wing_wet, horizontal_tail_wet, vertical_tail_wet, fuselage_wet)
            # print("the cd_0")
            # print(CD_0)
    
            # Use AVL to determine the CL_alpha of the wing based on the wing geometry -------------------------------------
            make_avl_file(c_root, c_tip, b, LE_sweep, dihedral, S, CD_0, inputs.M_cruise, avl_h, b_h, c_root_h, inputs.QC_sweep_h,
                          inputs.tap_h,
                          avl_v,
                          b_v, c_root_v, inputs.QC_sweep_v, inputs.tap_v, inputs.tail_type, 12, 5)
    
            new_Oswald, CD_cruise = run_avl(CL_cruise, inputs.M_cruise, CD_0)
    #        D = 0.5 * Rho_Cruise * (V_cruise ** 2) * S
            CL_alpha = (find_clalpha(inputs.M_cruise, CD_0, "conv_wing.avl") * 180) / np.pi
    
            iteration["updated parameters"] = [CD_0, inputs.Oswald, CL_alpha]
    
            # Determine maximum loads based on manoeuvring and gust envelope------------------------------------------------
            manoeuvring_loads = manoeuvring_envelope(W_TO, inputs.h_cruise, inputs.CL_Cruise_max, S, V_cruise)
            gust_loads = gust_envelope(W_TO, inputs.h_cruise, CL_alpha, S, mac, V_cruise, manoeuvring_loads[4])
    
            V_D = manoeuvring_loads[4][3]
    
            n_max_manoeuvring = max(manoeuvring_loads[1])
            n_max_gust = max(gust_loads[1])
    
            if n_max_manoeuvring > n_max_gust:
                n_ult = 1.5 * n_max_manoeuvring
    #            v_ult =1 manoeuvring_loads[0][list(manoeuvring_loads[1]).index(n_max_manoeuvring)]
            else:
                n_ult = 1.5 * n_max_gust
    #            v_ult = gust_loads[0][list(gust_loads[1]).index(n_max_gust)]
    
            # print(v_ult)
            # print(n_ult)
            # Note that for a conventional tail the horizontal tail starts at the root of the vertical tail, thus z_h=0
            if inputs.tail_type == 0:
                z_h = 0
            else:
                z_h = 2
    
            iteration["n_ult"] = n_ult
    
            # Start the class II weight estimation -------------------------------------------------------------------------
            # Determine the structural weight components -------------------------------------------------------------------
            # print("here")
            # print(t_over_c)
            # print(c_root)
            t_max_root = t_over_c * c_root
            
            weight_II = Class_II(W_TO, b, S, inputs.N_engines, V_D, T)
            w_weight = weight_II.wing_weight(W_F, HC_sweep, n_ult, t_max_root, inputs.wing_choice)
            mass_fractions[0] = (w_weight * lbs_to_kg * g_0) / W_TO
            # print("wing weight")
            # print(w_weight)
            emp_weight = weight_II.empennage_weight(inputs.empennage_choice, np.array([S_h, S_v]),
                                          np.array([HC_sweep_h, HC_sweep_v]),
                                          z_h, b_v)
            mass_fractions[1] = (emp_weight * lbs_to_kg * g_0) / W_TO
    
            # Note that the fuselage is assumed circular in this case
            fus_weight = weight_II.fuselage_weight(inputs.fuselage_choice, l_h, d_fuselage, d_fuselage, fuselage_wet)
            mass_fractions[2] = (fus_weight * lbs_to_kg * g_0) / W_TO
    
            # Note that the choice includes the engine choice here
            nac_weight = weight_II.nacelle_weight(inputs.nacelle_choice)
            mass_fractions[3] = (nac_weight * lbs_to_kg * g_0) / W_TO
    
            lg_weight = weight_II.landing_gear_weight()
            eng_weight = weight_II.engine_weight(inputs.w_engine)
    
            structural_weight = w_weight + emp_weight + fus_weight + nac_weight + lg_weight
    
            # Determine the propulsion system weight components ------------------------------------------------------------
            ai_weight = weight_II.induction_weight(inputs.duct_length, inputs.n_inlets, inputs.a_inlets, inputs.induction_choice)
    
            n_prop = inputs.prop_characteristics[0]
            n_blades = inputs.prop_characteristics[1]
            d_prop = inputs.prop_characteristics[2]
    
            if inputs.propeller_choice == 1:
                to_power = ((T * V_cruise) / (550 * inputs.prop_characteristics[3]))
                prop_weight = weight_II.propeller_weight(inputs.prop_choice, n_prop, d_prop, to_power, n_blades)
            else:
                prop_weight = 0
                to_power = 0
    
            # Note that choice is regarding type of fuel tanks
            fuel_sys_weight = weight_II.fuel_system_weight(inputs.n_fuel_tanks, W_F, inputs.fuel_sys_choice)
            # Choice is depending on type of engine controls and whether there is an afterburner
            w_ec = weight_II.calc_w_ec(l_fuselage, inputs.engine_choice)
            # Choice is depending on type of starting system and type of engine
            w_ess = weight_II.calc_w_ess(inputs.w_engine, inputs.start_up_choice)
            # Choice is depending on type of engine
            w_pc = weight_II.calc_w_pc(n_blades, n_prop, d_prop, to_power, inputs.prop_choice)
            # Choice is depending on type of engines
            w_osc = weight_II.calc_w_osc(inputs.oil_choice, inputs.w_engine)
    
            prop_sys_weight = eng_weight + ai_weight + prop_weight + fuel_sys_weight + w_ec + w_ess + w_pc + w_osc
            mass_fractions[4] = (prop_sys_weight * lbs_to_kg * g_0) / W_TO
    
            # Determine fixed equipment weight components ------------------------------------------------------------------
            # Calculate dynamic pressure at dive speed
            q_D = 0.5 * inputs.Rho_Cruise_input * V_D
            w_fc = weight_II.calc_w_fc(q_D)
    
            # Choice depends on type of engines
            w_hps_els = weight_II.calc_w_hps_els(W_E_I, inputs.hydro_choice, V_pax)
    
            # Maximum range has still to be determined
            w_instr = weight_II.calc_w_instr(W_E_I, inputs.maximum_range)
    
            w_api = weight_II.calc_w_api(V_pax, inputs.N_crew, inputs.N_pas)
            w_ox = weight_II.calc_w_ox(inputs.N_crew, inputs.N_pas)
            w_apu = weight_II.calc_w_apu()
            w_fur = weight_II.calc_w_fur(W_F)
            w_bc = weight_II.calc_w_bc(S_ff)
    
            fix_equip_weight = w_fc + w_hps_els + w_instr + w_api + w_ox + w_apu + w_fur + w_bc
            mass_fractions[5] = (fix_equip_weight * lbs_to_kg * g_0) / W_TO
    
            # Determine final operational empty weight ----------------------------
            W_E_II = (structural_weight + prop_sys_weight + fix_equip_weight) * lbs_to_kg * g_0
            
            # print(structural_weight)
            # print(prop_sys_weight)
            # print(fix_equip_weight)
            
            iteration["new empty weight"] = W_E_II
            W_e_frac = W_E_II / W_TO
    
            percentage = abs((W_E_II - W_E_I) / W_E_I)
    
            print("the percentage equals: " + str(percentage))
    #        print(iteration)
            percentages.append(percentage)
            empty_weight.append(W_E_II)
            i += 1
            
    return percentages
    
result = main_iterator(inputs.CL_cruise_input, inputs.CD_cruise_input, 
                       inputs.W_e_frac_input, inputs.fuel_fractions_input, 
                       inputs.mass_fractions_input, inputs.C_t0_input, inputs.V_cruise_input )
        

#CL_cruise, CD_cruise, W_e_frac, fuel_frac, mass_frac, C_t0, V_cruise
        
    # pie_chart_fracs.append(W_P / W_TO)
    # pie_chart_fracs.append(W_F / W_TO)
    # # pie_chart_fracs.append(W_E_II / W_TO)
    # pie_chart_fracs.append((fus_weight * lbs_to_kg * g_0) / W_TO)
    # pie_chart_fracs.append((emp_weight * lbs_to_kg * g_0) / W_TO)
    # pie_chart_fracs.append((w_weight * lbs_to_kg * g_0) / W_TO)
    # pie_chart_fracs.append((nac_weight * lbs_to_kg * g_0) / W_TO)
    # pie_chart_fracs.append((lg_weight * lbs_to_kg * g_0) / W_TO)
    # # pie_chart_fracs.append((eng_weight * lbs_to_kg * g_0) / W_TO)
    # pie_chart_fracs.append((prop_sys_weight * lbs_to_kg * g_0) / W_TO)
    # pie_chart_fracs.append((fix_equip_weight * lbs_to_kg * g_0) / W_TO)

    # print(fus_weight * lbs_to_kg)
    # print(emp_weight * lbs_to_kg)
    # print(prop_sys_weight*lbs_to_kg)

    # file = open("hbp_concept" + str(), "w")
    # file.write("The Mach number: " + str(M_cruise) + '\n')
    # file.write("The cruise altitude in m: " + str(h_cruise) + '\n')
    #
    # file.write("Take-off weight in N: " + str(round(W_TO, 2)) + '\n')
    # file.write("Empty weight in N: " + str(round(W_E_II, 2)) + '\n')
    # file.write("Wing weight in N: " + str(round(((w_weight * lbs_to_kg * g_0)), 2)) + '\n')
    # file.write("Fuselage weight in N: " + str(round(((fus_weight * lbs_to_kg * g_0)), 2)) + '\n')
    # file.write("Empennage weight in N: " + str(round(((emp_weight * lbs_to_kg * g_0)), 2)) + '\n')
    # file.write("Nacelle weight in N: " + str(round(((nac_weight * lbs_to_kg * g_0)), 2)) + '\n')
    # file.write("Engine weight in N: " + str(round(((eng_weight * lbs_to_kg * g_0)), 2)) + '\n')
    # file.write("Landing gear weight in N: " + str(round(((lg_weight * lbs_to_kg * g_0)), 2)) + '\n')
    # file.write("Fixed equipmentweight in N: " + str(round(((fix_equip_weight * lbs_to_kg * g_0)), 2)) + '\n')
    # file.write(
    #     "Propulsion system weight in N: " + str(round(((prop_sys_weight * lbs_to_kg * g_0)), 2)) + '\n')
    # file.write("Payload weight in N: " + str(round(W_P, 2)) + '\n')
    # file.write("Fuel weight in N: " + str(round(W_F, 2)) + '\n')
    # file.write("Thrust in N: " + str(round(T, 2)) + '\n')
    #
    # file.write("Cl_cruise: " + str(round(CL_cruise, 3)) + '\n')
    # file.write("Cd cruise: " + str(round(CD_cruise, 3)) + '\n')
    # file.write("Cruise density in kg/m^3: " + str(round(Rho_Cruise, 2)) + '\n')
    # file.write("Cd zero: " + str(round(CD_0, 4)) + '\n')
    # file.write("Lift over drag: " + str(round(CL_cruise / CD_cruise, 3)) + '\n')
    # file.write("Cruise speed in m/s: " + str(round(V_cruise, 2)) + '\n')
    #
    # file.write("The ultimate load factor: " + str(round(n_ult, 2)) + '\n')
    # file.write("The ultimate load factor speed in m/s: " + str(round(v_ult, 2)) + '\n')
    #
    # file.write("Span in m: " + str(round(b, 2)) + '\n')
    # file.write("Surface area in m^2: " + str(round(S, 2)) + '\n')
    # file.write("The length of the fuselage in m: " + str(l_fuselage) + '\n')
    # file.write("Root chord in m: " + str(round(c_root, 2)) + '\n')
    # file.write("The leading edge sweep in degrees: " + str(round(np.degrees(LE_sweep), 2)) + "\n")
    # file.write("The taper ratio: " + str(round(taper, 2)) + "\n")
    #
    # file.write("Aspect ratio: " + str(A) + '\n')
    # file.write("Oswald factor: " + str(round(Oswald, 2)) + '\n')
    # file.write("Wing weight in N: " + str(round(w_weight * lbs_to_kg * g_0, 2)) + '\n')
    # file.write("Number of engines: " + str(N_engines) + '\n')
    # file.write("Engine weight in N: " + str(w_engine * g_0) + "\n")
    #
    # file.write("Center of gravity locations in m: " + str(cg_locations) + "\n")
    # file.write("The weight fractions as % of take-off weight: " + str(pie_chart_fracs) + '\n')
    # file.close()

#        Velocity.append(V_cruise)
#        h.append(h_cruise)
#        AR.append(A)
#        Surface.append(S)
#        eff.append(Oswald)
#        cd_0.append(CD_0)
#        Ct0.append(Ct0_value)
#        Wcr.append(W_TO * 0.75)

# print(Velocity)
# print(h)
# print(AR)
# print(Surface)
# print(eff)
# print(cd_0)
# print(Ct0)
# print(Wcr)

#h_list, SAR_list = SAR(Velocity, h, AR, Surface, eff, cd_0, Ct0, Wcr)
#final_h.append(h_list)
#final_SAR.append(SAR_list)
#final_M.append(M_cruise)
#final_v.append(Velocity)
##
## pie_chart_labels = ["payload", "fuel", "fuselage", "empennage", "wing", "nacelle", "landing gear",
##                     "propulsion system", "fixed equipment"]
## print(len(pie_chart_labels))
## print(len(pie_chart_fracs))
## plt.pie(pie_chart_fracs, labels=pie_chart_labels, startangle=90, autopct='%.2f%%')
## plt.legend()
## plt.show()
#
#SAR_ref = 1.84236002771
#M_ref = 0.79
#H_ref = 12000
#min_SAR = []
#min_h = []
#h_v = []
#sar_v = []
##
#for r in range(len(final_h)):
#plt.plot(final_h[r], final_SAR[r], label='Mach %s' % round(final_M[r], 2))
#
## min_SAR.append(min(final_SAR[r]))
## i = final_SAR[r].index(min(final_SAR[r]))
## min_h.append(final_h[r][i])
#
## for q in range(len(final_v)):
##     for l in range(len(final_v[q])):
##         if final_v[q][l] > 700.:
##             ind = final_v[q].index(final_v[q][l])
##             h_v.append(final_h[q][ind])
##             sar_v.append(final_SAR[q][ind])
##             break
#
## plt.plot(h_v, sar_v, label="Minimum required speed for cost")
#plt.plot(H_ref, SAR_ref / 200, 'mo', label="Ref. aircraft")
## plt.plot(min_h, min_SAR, label="Optimum line")
#plt.hlines(0.9 * SAR_ref / 200, 0, 13000, "gray", "--")
#plt.legend()
#plt.title('Fuel consumption per passenger w.r.t. Mach number')
#plt.xlabel("Altitude [m]")
## plt.ylim(0.006, 0.011)
#plt.ylabel("Fuel consumption [kg/km/passenger]")
## plt.savefig("Design 1")
#plt.show()

# final_diagram(CD_0, Oswald)
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

# iterations = np.arange(0, len(percentages), 1)
# # final_diagram(CD_0, Oswald)
# plt.figure(1)
# plt.plot(iterations, percentages)
# plt.figure(2)
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
