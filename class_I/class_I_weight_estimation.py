from constants_and_conversions import g_0, per_hr_to_N
import numpy as np
#import matplotlib.pyplot as plt


#   calculate the total payload weight
def calc_payload_weight(n_passengers, n_crew, m_person, m_cargo):
    w_people = (n_passengers + n_crew) * m_person * g_0
    w_cargo = n_passengers * m_cargo * g_0
    
    return w_people + w_cargo


#   calculate the cruise coefficent based on ADSEE
def calc_cruise_coefficient(cl, cd, r, V, c_j):
#    print("The range equals " + str(r))
#    print("The V equals "+ str(V))
#    print("The lift over drag ratio equals " + str(cl/cd))
#    print(r, cl, cd, V, c_j)
    ff = 1 / (np.e ** (r / ((V / (c_j * g_0)) * (cl / cd))))
#    print(ff)
    return ff
#    return 1 / (np.e ** (r / ((V / (c_j * per_hr_to_N * g_0)) * (cl / cd))))


#   calculate the loiter coefficient based on ADSEE
def calc_loiter_coefficient(A, oswald_factor, cd_0, E, c_j):
    cd = 2 * cd_0
    cl = np.sqrt((np.pi * A * oswald_factor * cd_0))
    ff = 1 / (np.e ** (E / ((1 / (c_j * g_0)) * (cl / cd))))
    return ff
    #return 1 / (np.e ** (E / ((1 / (c_j * g_0 * per_hr_to_N)) * (cl / cd))))


#   calculate the total required fuel by multiplying all phases of the flight profile
def calc_fuel_fraction(coefficients):
    fraction = 1
    
    for number in coefficients:
        fraction = fraction * number
        
    return 1 - fraction


#   final calculation of class I weight estimation
def class_I(cl, cd, r_cruise, r_res, v_cruise, cj_cruise, W_tfo_frac, W_e_frac, 
            fuel_fractions, N_pas, N_crew, M_person, M_cargo):  
    
    cruise_1 = calc_cruise_coefficient(cl, cd, r_cruise, v_cruise, cj_cruise)
#    loiter = calc_loiter_coefficient(A, Oswald, S_ratio, C_fe, E, c_j_loiter)
    cruise_2 = calc_cruise_coefficient(cl, cd, r_res, v_cruise, cj_cruise)

#    print("The first cruise coefficient equals " + str(cruise_1))
#    print("The loiter coefficient equals " + str(loiter))
#    print("The second cruise coefficient equals " + str(cruise_2))

    fuel_fractions[-1] = cruise_1
    fuel_fractions[-2] = cruise_2
#    fuel_fractions = np.append(fuel_fractions, [cruise_1, cruise_2])
#    print(fuel_fractions)
    W_f_frac = calc_fuel_fraction(fuel_fractions)
#    fuel_fractions = [start, taxi, t_o, climb_1, descent_1, climb_2, 
#                        descent_2, landing, 0, 0]
    nom_fuel_fractions = fuel_fractions[0:5] + [fuel_fractions[7]] + [cruise_1]
    W_nom_fuel_frac = calc_fuel_fraction(nom_fuel_fractions)
    
    W_to_frac = 1 - W_f_frac - W_tfo_frac - W_e_frac
    
#    print(W_f_frac) 
#    print(W_e_frac)
#    print(W_f_frac)
    W_P = calc_payload_weight(N_pas, N_crew, M_person, M_cargo)
    W_TO = W_P / W_to_frac
    W_F = W_f_frac * W_TO
    W_E = W_e_frac * W_TO + W_tfo_frac * W_TO
    W_nom_F = W_nom_fuel_frac * W_TO
    return np.array([W_TO, W_E, W_P, W_F]), W_nom_F

# r_cruise = np.array([1200000, 1400000, 1600000, 1800000])
# results = np.zeros((4, len(r_cruise)))
# percentages = np.zeros((3, len(r_cruise) - 1))
# for i in range(len(r_cruise)):
#     result = class_I(CL_cruise_input, CD_cruise_input, r_cruise[i], reserve_range, V_cruise, c_j_cruise, W_tfo_frac, W_e_frac_input,
#                      fuel_fractions, N_pas, N_crew, W_person)
#     for j in range(4):
#         results[j][i] = result[j]
#
#     if i > 0:
#         percentages[0][i - 1] = ((r_cruise[i] - r_cruise[i - 1]) / r_cruise[i - 1]) * 100
#         percentages[1][i - 1] = ((result[1] - results[1][i - 1]) / results[1][i - 1]) * 100
#         percentages[2][i - 1] = ((result[3] - results[3][i - 1]) / results[3][i - 1]) * 100
#
# # plt.plot(percentages[0], percentages[1])
# # plt.plot(r_cruise, results[3])
# # plt.plot(r_cruise, results[2])
# # plt.plot(r_cruise, results[1])
# # plt.show()
#
# print(results[3], r_cruise)
# narrow_body = class_I(15, 1, 3550 * nm_to_km * 1000, 250 * nm_to_km * 1000, 450 * kts_to_ms, 0.4, 0.003, 0.48,
#                       fuel_fractions_input, 180, 5, W_person, W_carg, 20882*9.80565)
# wide_body = class_I(17, 1, 4820 * nm_to_km * 1000, 250 * nm_to_km * 1000, 476 * kts_to_ms, 0.25, 0.003, 0.56,
#                       fuel_fractions_input, 440, 11, W_person, W_carg, 54635*9.80565)
#
# print(narrow_body)
# print(wide_body)
