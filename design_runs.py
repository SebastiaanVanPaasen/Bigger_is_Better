from read_csv_input import read_input
from iterator import main_iterator

import numpy as np

def run_designs(design_names):
    W_TO, W_P, W_F, W_E_II, w_weight, emp_weight, fus_weight, nac_weight, prop_weight, fix_equip_weight = ["Take off weight [N]"], ["Payload weight [N]"], ["Fuel weight [N]"], ["Empty weight [N]"], ["Wing weight [N]"], ["Empennage weight [N]"], ["Fuselage weight [N]"], ["Nacelle weight [N]"], ["Propulsion system weight [N]"], ["Fixed equipment weight [N]"]
    weight_results = [W_TO, W_P, W_F, W_E_II, w_weight, emp_weight, fus_weight, nac_weight, prop_weight, fix_equip_weight]
    
    CL_cr, CD_cr, LD, c_t, Oswald = ["CL cruise [-]"], ["CD cruise [-]"], ["L/D-ratio"], ["C_t [kg/Ns]"], ["Oswald factor [-]"]
    coefficient_results = [CL_cr, CD_cr, LD, c_t, Oswald]
    
    A, S, b, c_root, c_tip, QC_sweep, taper = ["Aspect ratio [-]"], ["Surface area [m^2]"], ["Span [m]"], ["Root chord [m]"], ["Tip chord [m]"], ["Quarter chord sweep [deg]"], ["Taper ratio [-]"]
    planform_results = [A, S, b, c_root, c_tip, QC_sweep, taper]
    
    l_fuselage, d_fuselage, l_nosecone, l_tailcone, l_h = ["Fuselage length [m]"], ["Fuselage diameter [m]"], ["Nosecone length [m]"], ["Tailcone length [m]"], ["Tail arm [m]"]
    fuselage_results = [l_fuselage, d_fuselage, l_nosecone, l_tailcone, l_h]
    
    S_h, S_v = ["Horizontal tail surface area [m^2]"], ["Vertical surface area [m^2]"]
    tail_results = [S_h, S_v]
    
    T_TO, M_cr, V_cr, H_cr, sar = ["Take off thrust [N]"], ["Cruise Mach number [-]"], ["Cruise velocity [m/s]"], ["Cruise altitude [m]"], ["Specific fuel consumption [kg/km/pas]"]
    env_results = [T_TO, M_cr, V_cr, H_cr, sar]
    
    F_decr, RF_decr, C_decr = ["Fuel consumption reduction [%]"], ["RF reduction [%]"], ["DOC reduction [%]"]
    diff_results = [F_decr, RF_decr, C_decr]
    
    names = [""] + design_names
    finals = []
    
    for i in range(len(design_names)): 
        print("Starting on design: " + str(design_names[i]))
        cf, char, env, eng, opti, tail = read_input(design_names[i])
        weights, coefficients, planform, fus, tails, environment, diff = main_iterator(cf, char, env, eng, opti, tail)
    
        for i in range(len(weights)):
            weight_results[i].append(weights[i])
        finals.append(weight_results)
            
        for i in range(len(coefficients)) :
            coefficient_results[i].append(coefficients[i])
        finals.append(coefficient_results)
        
        for i in range(len(planform_results)):
            planform_results[i].append(planform[i])
        finals.append(planform_results)
            
        for i in range(len(fuselage_results)):
            fuselage_results[i].append(fus[i])
        finals.append(fuselage_results)
            
        for i in range(len(tail_results)):
            tail_results[i].append(tails[i])
        finals.append(tail_results)
            
        for i in range(len(env_results)):
            env_results[i].append(environment[i])
        finals.append(env_results)
        
        for i in range(len(diff)):
            diff_results[i].append(diff[i])
        finals.append(diff_results)
        
    final_result = np.asarray(names)
    
#    print(final_result)
    for i in range(len(finals)):
        final_result = np.append(final_result, np.asarray(finals[i]))
        
    final_result = np.reshape(final_result, (-1, (len(design_names)+1)))
    
    np.savetxt("csv_results_LOW_2E_SEMIDD2.csv", final_result, fmt= '%s', delimiter=";")

designs = []
for i in range(54):
    designs.append("Design " + str(i+1))
    
#print(designs)
    
run_designs(designs)
