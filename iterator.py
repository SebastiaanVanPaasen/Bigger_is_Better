import input_files.standard_inputs as ip
import numpy as np

from constants_and_conversions import Temp_0, g_0, gamma, a, R_gas, Rho_0, lbs_to_kg
from class_II_weight_estimation import Class_II
from class_I.class_I_weight_estimation import class_I
from class_I.fuselage_cross_section import fuselage_cross_section
from class_I.planform import wing_parameters, determine_half_chord_sweep
from class_I.drag_estimation import Wing_wetted_area, H_tail_wetted_area, V_tail_wetted_area, Fus_wetted_area, \
    Zero_Lift_Drag_est
from class_I.class_I_empennage_landinggear import class_I_empennage  
#, _calc_h_tail_II, _calc_v_tail_II
from class_I.flight_envelope import manoeuvring_envelope, gust_envelope
from avl.conv_wing_avl import make_avl_file, run_avl, find_clalpha
from RF_calc import Radiative, Costs
#from class_I.seats import cg_seats, W_seats
#from class_I.loading_diagram import potato
#from sc_sensitivity import class_II_empennage
#from stability_and_control.stability import C_L_alpha_Ah, Sh_S_stability
#from stability_and_control.control_curve import Sh_S_control
#from stability_and_control.control_stability import control_stability_plot


def main_iterator(cf, char, env, eng, opt, tails):
    CL_cr, CD_cr, c_t0, Oswald = cf[0], cf[1], cf[2], cf[3]
    T_TO_ip, S_ip, A, wing, tail = char[0], char[1], char[2], char[3], char[4]
    V_ref, H_cr, M_cr = env[0], env[1], env[2]
    N_eng, m_eng, l_inl, A_inl = eng[0], eng[1], eng[2], eng[3]
    A_h, QC_sweep_h, tap_h = tails[0], tails[1], tails[2]
    A_v, QC_sweep_v, tap_v = tails[3], tails[4], tails[5]
    N_pas_below = tails[6]
    
    if H_cr < 11000.:
        Temp_cr = Temp_0 + a * H_cr  # K  based on the altitude you fly at
    else:
        Temp_cr = 216.65
    
    a_cr = np.sqrt(gamma * R_gas * Temp_cr)  # m/s based on the temperature
    
    Rho_cr = Rho_0 * ((1 + (a * H_cr) / Temp_0) ** (-(g_0 / (R_gas * a))))  # kg/m^3   based on cruise altitude
    
    V_cr = M_cr * a_cr
    
    c_t = (c_t0 / V_ref) * V_cr
    
#    if N_eng == 2:
#        b_n_1, b_n_2, b_n_3, b_n_4 = 3.2, 3.2, 0, 0
#        l_n_1, l_n_2, l_n_3, l_n_4 = -2, -2, 0, 0
#        
#    else:
#        b_n_1, b_n_2, b_n_3, b_n_4 = 3.2, 3.2, 3.2, 3.2
#        l_n_1, l_n_2, l_n_3, l_n_4 = -2, -2, -2, -2
    
    mass_fractions = ip.mass_fractions
    percentages = []

    fuselage_design = fuselage_cross_section(ip.N_pas, N_pas_below)
    
    # Define fuselage parameters out of the fuselage design
    d_fuselage = fuselage_design[1]
    l_nosecone = np.sum(fuselage_design[6]) / 2
    l_cabin = fuselage_design[2]
    l_tailcone = np.sum(fuselage_design[5]) / 2
    l_fuselage = fuselage_design[7]
    V_pax = 0.25 * np.pi * (d_fuselage ** 2) * l_cabin  
    # m^3  Volume of the passenger cabin

#    xcg_seats = cg_seats(fuselage_design[16], fuselage_design[11], fuselage_design[12],l_nosecone)
    
#    W_window, W_aisle, W_middle = W_seats(fuselage_design[8], fuselage_design[11], fuselage_design[12], fuselage_design[15])


#    emp_constants = [ip.N_cargo, l_fuselage, ip.cargo_fwdfrac, xcg_seats, W_window, W_aisle, W_middle, ip.N_pas, g_0*ip.m_person]
#    
#    M_h_cr = M_cr * ip.V_h_norm
#    M_w_landing = ip.V_stall_l / np.sqrt(gamma * R_gas * Temp_cr)
    
    # Define fineness ratio's based on the fuselage design
#    f = l_fuselage / d_fuselage
#    f_tc = l_tailcone / d_fuselage
#    f_nc = l_nosecone / d_fuselage

    # Area of freight floor, very rough initial guess for now
    S_ff = 0.5 * d_fuselage * l_cabin
    
    # Starting the iteration process --------------------------------------
    i = 0
    maximum = 20
    percentage = 1
    W_e_frac = ip.w_e
    while i < maximum and percentage > 0.01:
        print("Starting on iteration: " + str(i))
        # Performing class I weight estimation ----------------------------
        weights = class_I(CL_cr, CD_cr, ip.m_range, ip.r_range, V_cr, c_t, 
                          ip.w_tfo, W_e_frac, ip.fuel_fractions, ip.N_pas, 
                          ip.N_crew, ip.m_person, ip.m_carg)

        W_TO, W_E_I, W_P, W_F = weights[0], weights[1], weights[2], weights[3]
        
        # print(W_TO)
        # print(V_cruise)
        # print(W_tfo_frac)
        
        mass_fractions[7] = (weights[1] / weights[0])  # empty mass fraction
        mass_fractions[8] = (weights[2] / weights[0])  # payload mass fraction
        mass_fractions[9] = (weights[3] / weights[0])  # fuel mass fraction

        # Choosing a design point based on the T/W-W/S diagram ------------
        T_TO, S = W_TO * T_TO_ip, W_TO / S_ip
        CL_cr = (0.75 * W_TO) / (0.5 * Rho_cr * (V_cr ** 2) * S)
        if CL_cr > 1.:
            CL_cr = 1.
            print("The CL is 1")

#            print(Rho_Cruise, V_cruise, S, CD_cruise)

        # Sizing the wing based on cruise parameters and wing configuration 
        QC_sweep, LE_sweep, b, c_root, c_tip, di, t_c, mac, taper = wing_parameters(M_cr, CL_cr, S, A, wing)

        # Calculate additional required half chord sweep for later use
        HC_sweep = determine_half_chord_sweep(c_tip, QC_sweep, c_root, b)


        # Perform first order cg-range estimation based on statistics -----
        x_payload = 0.5 * l_fuselage  # m     cg-location payload w.r.t. nose
        
#        if i == 0:
        cg_locations, tail_h, tail_v, x_lemac, avl_h, avl_v, xcg_aft = class_I_empennage(mass_fractions, mac, l_fuselage,
                                                                            ip.xcg_eng, ip.l_nac, 
                                                                            ip.xcg_oew_mac, x_payload,
                                                                            d_fuselage, b, S, taper, 
                                                                            LE_sweep, [A_v, tap_v, QC_sweep_v], 
                                                                            [A_h, tap_h, QC_sweep_h])
         
        l_h, c_root_h, c_tip_h, b_h, S_h = tail_h[0], tail_h[1], tail_h[2], tail_h[3], tail_h[4]
        c_root_v, c_tip_v, b_v, S_v = tail_v[0], tail_v[1], tail_v[2], tail_v[3] 
#            opt_Sh_S = S_h / S
#            opt_X_LEMAC = x_lemac 
                              
            
#        else:
#            
#            c_f = 0.25 * mac #flap chord
#            b_f0 = 0.4 * b #end of flap span
#            b_fi = 0.1 * b #begin of flap span
#            
##            old_S_h = S_h
#            S_h = opt_Sh_S * S
##            S_v = S_h / old_S_h * S_v
#            
#            x_lemac = opt_X_LEMAC
##            print("The leading edge of mac " + str(x_lemac))
#            tail_h = _calc_h_tail_II(xcg_aft, A_h, l_fuselage, tap_h, QC_sweep_h, S_h)
#            tail_v = _calc_v_tail_II(A_v, l_fuselage, tap_v, S_v, QC_sweep_v)
#            min_cg, max_cg, X_LEMAC_range, min_cg_range = potato(l_nosecone, W_TO, ip.xcg_eng, ip.l_nac, mass_fractions, tail_h, 
#                                                         tail_v, ip.s_m, cg_locations, x_lemac, emp_constants, mac)
#            
#            Y_MAC = (b / 6.) * ((1 + 2 * taper) / (1 + taper))
#            
#            x_cg = np.linspace(min(min_cg), max(max_cg))
#            x_le_h, sweep_LE_h, y_MAC_h, MAC_h, l_h = tail_h
#            x_le_v, sweep_LE_v, y_MAC_v, MAC_v = tail_v
#            HC_sweep_h = determine_half_chord_sweep(c_tip_h, QC_sweep_h, c_root_h, b_h)
#            x_LE_root = x_lemac - LE_sweep * Y_MAC
#            
#            if wing == 0: 
#                if tail == 0:
#                    
##                    print("vertical tail part " + str(np.tan(sweep_LE_v) * b_v))
#                    m_tv = fuselage_design[1] / 2.
##                    print("the new tail arm " + str(l_h))
#    
#                else:
#                    m_tv = fuselage_design[1] / 2. + b_v
#                    
##                    print("vertical tail part " + str(np.tan(sweep_LE_v) * b_v))
#                    l_h = l_h + (x_le_h - x_le_v) + np.tan(sweep_LE_v) * b_v
##                    print("The new tail arm " + str(l_h))
#                    
#                    
#                    
#                
#            else:
#                if tail == 0:
#                    m_tv = - fuselage_design[1] / 2.
##                    print("the new tail arm " + str(l_h))
##                    print("vertical tail part " + str(np.tan(sweep_LE_v) * b_v))
#                else:
##                    print("the previous tail arm " + str(l_h))
##                    print("vertical tail part " + str(np.tan(sweep_LE_v) * b_v))
#                    m_tv = b_v - fuselage_design[1] / 2.
#                    l_h = l_h + (x_le_h - x_le_v) + np.tan(sweep_LE_v) * b_v
##                    print("The new tail arm " + str(l_h))
#
#                    
#            
#            stability_lessmargin_list, stability_list = Sh_S_stability(x_cg, M_h_cr, ip.eta, HC_sweep_h, HC_sweep, A, A_h, M_cr, fuselage_design[1], b, S,
#                            l_h, QC_sweep, m_tv, ip.V_h_norm, b_n_1, b_n_2, b_n_3, b_n_4, l_n_1, l_n_2, l_n_3, l_n_4, ip.x_ac_wing, 
#                            fuselage_design[1], x_LE_root, mac, c_root, mac, taper, ip.SM)
#            
#              
#        
#            
#            
#            C_L_alpha_Ah_landing = C_L_alpha_Ah(M_w_landing, ip.eta, HC_sweep, A, fuselage_design[1], b, S, c_root)  
#            
#            
#            C_L_Ah_landing = C_L_alpha_Ah_landing * ip.alpha_land
#            
#            
#            
#            S_netfus = S - (fuselage_design[1] * c_root)
#            control_list = Sh_S_control(ip.CL_H, C_L_Ah_landing, l_h, ip.V_h_norm, x_cg, ip.Cm_0, A, QC_sweep,
#                         ip.delta_flap ,b ,b_f0 ,b_fi ,taper , mac, c_f, ip.dc_c_f, ip.mu_2, ip.mu_3, ip.x_ac, ip.CL_l_max,
#                         fuselage_design[1], fuselage_design[1], fuselage_design[7],S,S_netfus,HC_sweep, M_w_landing, ip.eta, ip.CL_0)
#            
#            
#            
#            opt_X_LEMAC, opt_Sh_S, xcg_aft = control_stability_plot(x_cg, min_cg, max_cg, X_LEMAC_range, control_list, stability_list, mac, l_fuselage)

#        print("tail surface area " + str(S_h/S))
        # Calculate the accompanying tail sizes ---------------------------
        HC_sweep_h = determine_half_chord_sweep(c_tip_h, QC_sweep_h, c_root_h, b_h)
        HC_sweep_v = determine_half_chord_sweep(c_tip_v, QC_sweep_v, c_root_v, b_v)

        # Calculate new wetted area----------------------------------------
        main_wing_wet = Wing_wetted_area(c_root, c_tip, d_fuselage, b, S)
        horizontal_tail_wet = H_tail_wetted_area(c_root_h, tap_h, b_h)
        vertical_tail_wet = V_tail_wetted_area(c_root_v, tap_v, b_v)
        fuselage_wet = Fus_wetted_area(d_fuselage, l_nosecone, l_cabin, l_tailcone)

        # Determine the new CD_0 value ------------------------------------
        CD_0 = Zero_Lift_Drag_est(S, main_wing_wet, horizontal_tail_wet, vertical_tail_wet, fuselage_wet)
        
#            print("the cd_0")
#            print(CD_0)

        # Use AVL to determine the CL_alpha of the wing based on the wing geometry -------------------------------------
        make_avl_file(c_root, c_tip, b, LE_sweep, di, S, CD_0, M_cr, avl_h,
                      b_h, c_root_h, QC_sweep_h, tap_h, avl_v, 
                      b_v, c_root_v, QC_sweep_v, tap_v, tail, 12, 5)

        new_Oswald, CD_cr = run_avl(CL_cr, M_cr, CD_0)
        CL_alpha = (find_clalpha(M_cr, CD_0, "conv_wing.avl") * 180) / np.pi

        # Determine maximum loads based on manoeuvring and gust envelope------------------------------------------------
        manoeuvring_loads = manoeuvring_envelope(W_TO, H_cr, ip.CL_cr_max, S, V_cr)
        gust_loads = gust_envelope(W_TO, H_cr, CL_alpha, S, mac, V_cr, manoeuvring_loads[4])

        V_D = manoeuvring_loads[4][3]

        n_max_manoeuvring = max(manoeuvring_loads[1])
        n_max_gust = max(gust_loads[1])

        if n_max_manoeuvring > n_max_gust:
            n_ult = 1.5 * n_max_manoeuvring
#            v_ult =1 manoeuvring_loads[0][list(manoeuvring_loads[1]).index(n_max_manoeuvring)]
        else:
            n_ult = 1.5 * n_max_gust
#            v_ult = gust_loads[0][list(gust_loads[1]).index(n_max_gust)]

#            print(v_ult)
#            print(n_ult)

        # Note that for a conventional tail the horizontal tail starts at 
        # the root of the vertical tail, thus z_h=0
        if tail == 0:
            z_h = 0
        else:
            z_h = 2

        # Start the class II weight estimation ----------------------------
        # Determine the structural weight components ----------------------
        t_max_root = t_c * c_root
        
        weight_II = Class_II(W_TO, b, S, N_eng, m_eng, V_D, T_TO)
        w_weight = weight_II.wing_weight(W_F, HC_sweep, n_ult, t_max_root, opt[0]) * lbs_to_kg * g_0
        mass_fractions[0] = w_weight / W_TO

        emp_weight = weight_II.empennage_weight(opt[1], np.array([S_h, S_v]),
                                      np.array([HC_sweep_h, HC_sweep_v]),
                                      z_h, b_v) * lbs_to_kg * g_0
        mass_fractions[1] = emp_weight / W_TO

        # Note that the fuselage is assumed circular in this case
        fus_weight = weight_II.fuselage_weight(opt[2], l_h, d_fuselage, d_fuselage, fuselage_wet) * lbs_to_kg * g_0
        mass_fractions[2] = fus_weight / W_TO

        # Note that the choice includes the engine choice here
        nac_weight = weight_II.nacelle_weight(opt[3]) * lbs_to_kg * g_0
        mass_fractions[3] = nac_weight / W_TO

        lg_weight = weight_II.landing_gear_weight() * lbs_to_kg * g_0
        mass_fractions[6] = lg_weight / W_TO

        structural_weight = w_weight + emp_weight + fus_weight + nac_weight + lg_weight

        # Determine the propulsion system weight components ---------------
        eng_weight = weight_II.engine_weight() * lbs_to_kg * g_0
                    
        ai_weight = weight_II.induction_weight(l_inl, N_eng, A_inl, opt[4]) * lbs_to_kg * g_0


        if opt[5] == 1:
            n_prop = 4
            n_blades = 800
            d_prop = 3
            to_power = ((T_TO* V_cr) / (550 * 10))
            prop_weight = weight_II.propeller_weight(n_prop, d_prop, to_power, n_blades) * lbs_to_kg * g_0
            w_pc = weight_II.calc_w_pc(n_blades, n_prop, d_prop, to_power, opt[5]) * lbs_to_kg * g_0
        else:
            prop_weight = 0
            to_power = 0
            w_pc = 0

        # Note that choice is regarding type of fuel tanks
        fuel_sys_weight = weight_II.fuel_system_weight(ip.N_tanks, W_F, opt[6]) * lbs_to_kg * g_0
        # Choice is depending on type of engine controls and whether there is an afterburner
        w_ec = weight_II.calc_w_ec(l_fuselage, opt[6]) * lbs_to_kg * g_0
        # Choice is depending on type of starting system and type of engine
        w_ess = weight_II.calc_w_ess(opt[7][0]) * lbs_to_kg * g_0
         # Choice is depending on type of engines
        w_osc = weight_II.calc_w_osc(opt[8]) * lbs_to_kg * g_0

        prop_sys_weight = eng_weight + ai_weight + prop_weight + fuel_sys_weight + w_ec + w_ess + w_pc + w_osc
        mass_fractions[4] = prop_sys_weight / W_TO

        # Determine fixed equipment weight components ---------------------        
        # Calculate dynamic pressure at dive speed
        q_D = 0.5 * Rho_cr * V_D
        w_fc = weight_II.calc_w_fc(q_D) * lbs_to_kg * g_0

        # Choice depends on type of engines
        w_hps_els = weight_II.calc_w_hps_els(W_E_I, opt[9], V_pax) * lbs_to_kg * g_0

        # Maximum range has still to be determined
        w_instr = weight_II.calc_w_instr(W_E_I, ip.max_range) * lbs_to_kg * g_0

        w_api = weight_II.calc_w_api(V_pax, ip.N_crew, ip.N_pas) * lbs_to_kg * g_0
        w_ox = weight_II.calc_w_ox(ip.N_crew, ip.N_pas) * lbs_to_kg * g_0
        w_apu = weight_II.calc_w_apu() * lbs_to_kg * g_0
        w_fur = weight_II.calc_w_fur(W_F) * lbs_to_kg * g_0
        w_bc = weight_II.calc_w_bc(S_ff) * lbs_to_kg * g_0

        fix_equip_weight = w_fc + w_hps_els + w_instr + w_api + w_ox + w_apu + w_fur + w_bc
        mass_fractions[5] = fix_equip_weight / W_TO

        # Determine final operational empty weight ----------------------------
        W_E_II = structural_weight + prop_sys_weight + fix_equip_weight
        
        W_TO = W_E_II + W_P + W_F
        W_e_frac = W_E_II / W_TO

        percentage = abs((W_E_II - W_E_I) / W_E_I)

#        print("the percentage equals: " + str(percentage))
#        print(iteration)
        percentages.append(percentage)
        
        i += 1
    
    sar = ((((0.5 * Rho_cr * (V_cr ** 2) * S * CD_cr) * c_t) / V_cr) * 1000 )/ ip.N_pas  
    # in kg/km/pas
#    f_decr, rf_decr, cost_decr 
    RF_decr, F_decr = Radiative(H_cr, sar)
    C_decr = Costs(V_cr, sar)
#    print(values) 
#    print(F_decr)
#    print(C_decr)
#    print(RF_decr)
#    print(f_decr, rf_decr, cost_decr)
    
    weights = [W_TO, W_P, W_F, W_E_II, w_weight, emp_weight, fus_weight, nac_weight, prop_sys_weight, fix_equip_weight]
    coefficients = [CL_cr, CD_cr, CL_cr/CD_cr, c_t, Oswald]
    planform = [A, S, b, c_root, c_tip, QC_sweep, taper]
    fus = [l_fuselage, d_fuselage, l_nosecone, l_tailcone, l_h]
    tails = [S_h, S_v]
    environment = [T_TO, M_cr, V_cr, H_cr, sar]
    differences = [F_decr, RF_decr, C_decr]
    
    return weights, coefficients, planform, fus, tails, environment, differences


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

