from read_csv_input import read_input
from iterator import main_iterator
import matplotlib.pyplot as plt
import numpy as np

def run_designs(design_names):
    W_TO, W_P, W_F, W_E_II, w_weight, emp_weight, fus_weight, nac_weight, prop_weight, fix_equip_weight = ["Take off weight [N]"], ["Payload weight [N]"], ["Fuel weight [N]"], ["Empty weight [N]"], ["Wing weight [N]"], ["Empennage weight [N]"], ["Fuselage weight [N]"], ["Nacelle weight [N]"], ["Propulsion system weight [N]"], ["Fixed equipment weight [N]"]
    weight_results = [W_TO, W_P, W_F, W_E_II, w_weight, emp_weight, fus_weight, nac_weight, prop_weight, fix_equip_weight]
    
    CL_cr, CD_cr, LD, c_t, Oswald, CD_0 = ["CL cruise [-]"], ["CD cruise [-]"], ["L/D-ratio"], ["C_t [kg/Ns]"], ["Oswald factor [-]"], ["CD_0 [-]"]
    coefficient_results = [CL_cr, CD_cr, LD, c_t, Oswald, CD_0]
    
    A, S, b, c_root, c_tip, QC_sweep, taper, mac, di, t_c, x_LE_root = ["Aspect ratio [-]"], ["Surface area [m^2]"], ["Span [m]"], ["Root chord [m]"], ["Tip chord [m]"], ["Quarter chord sweep [deg]"], ["Taper ratio [-]"], ["Mean aerodynamic chord [m]"], ["Dihedral [deg]"], ["Thickness over chord [-]"], ["Leading edge root location [m]"]
    planform_results = [A, S, b, c_root, c_tip, QC_sweep, taper, mac, di, t_c, x_LE_root]
    
    l_fuselage, d_fuselage, l_nosecone, l_tailcone, l_cabin, _min_, _max_ = ["Fuselage length [m]"], ["Fuselage outer diameter [m]"], ["Nosecone length [m]"], ["Tailcone length [m]"], ["Cabin length [m]"], ["Forward cg [%MAC]"], ["Aft cg [%MAC]"]
    fuselage_results = [l_fuselage, d_fuselage, l_nosecone, l_tailcone, l_cabin, _min_, _max_]
    
    S_h, x_le_h, sweep_LE_h, MAC_h, l_h, c_root_h, c_tip_h, b_h, S_v, x_le_v, sweep_LE_v, MAC_v, c_root_v, c_tip_v, b_v = ["Horizontal tail surface area [m^2]"], ["Horizontal tail leading edge root position [m]"], ["Horizontal tail leading edge sweep [deg]"], ["Horizontal tail mean aerodynamic chord [m]"], ["Horizontal tail arm [m]"], ["Horizontal tail root chord [m]"], ["Horizontal tail tip chord [m]"], ["Horizontal tail span [m]"], ["Vertical tail surface area [m^2]"], ["Vertical tail leading edge root position [m]"], ["Vertical tail leading edge sweep [deg]"], ["Vertical tail mean aerodynamic chord [m]"], ["Vertical tail root chord [m]"], ["Vertical tail tip chord [m]"], ["Vertical tail span [m]"]
    tail_results = [S_h, x_le_h, sweep_LE_h, MAC_h, l_h, c_root_h, c_tip_h, b_h, S_v, x_le_v, sweep_LE_v, MAC_v, c_root_v, c_tip_v, b_v]
    
    T_TO, M_cr, V_cr, H_cr, sar, Tot_fuel_cons = ["Take off thrust [N]"], ["Cruise Mach number [-]"], ["Cruise velocity [m/s]"], ["Cruise altitude [m]"], ["Cruise fuel consumption [kg/km/pas]"], ["Mission fuel consumption [kg/km/pas]"]
    env_results = [T_TO, M_cr, V_cr, H_cr, sar, Tot_fuel_cons]
    
    F_decr, Tot_fuel_diff, RF_decr, C_decr = ["Fuel consumption reduction [%]"], ["Total fuel consumption reduction [%]"], ["RF reduction [%]"], ["DOC reduction [%]"]
    diff_results = [F_decr, Tot_fuel_diff, RF_decr, C_decr]
    
#    X_LE_ROOT, l_cabin_above, wing_gap = ["Wing position root [m]"], ["Length upper deck [m]"], ["Gap between wing and cabin [m]"]
#    wingpos_results = [X_LE_ROOT, l_cabin_above, wing_gap]
    
    names = [""] + design_names
    finals = []
    above_list = []
    below_list = []
    pas_list = []
    
    for i in range(len(design_names)): 
        print("Starting on design: " + str(design_names[i]))
        cf, char, env, eng, opti, tail = read_input(design_names[i])
        weights, coefficients, planform, fus, tails, environment, diff = main_iterator(cf, char, env, eng, opti, tail)
        
        above_list.append(fus[-3])
        below_list.append(fus[-2])
        pas_list.append(fus[-1])
        
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
        
#        for i in range(len(wingpos)):
#            wingpos_results[i].append(wingpos[i])
#        finals.append(wingpos_results)
        
    final_result = np.asarray(names)
#    plt.plot(pas_list, above_list, 'rx')
#    plt.plot(pas_list, below_list, 'bx')
#    plt.show()

#    print(final_result)
    for i in range(len(finals)):
        final_result = np.append(final_result, np.asarray(finals[i]))
        
    final_result = np.reshape(final_result, (-1, (len(design_names)+1)))
    
    np.savetxt("test1.csv", final_result, fmt= '%s', delimiter=";")

designs = []
for i in range(1):
    designs.append("Design " + str(i+1))
    
#print(designs)
      
run_designs(designs)     

