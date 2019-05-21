from constants_and_conversions import *
import numpy as np
# from input_files.strutted_wing import *
import matplotlib.pyplot as plt


def calc_payload_weight(n_passengers, n_crew, w_person, w_cargo):
    # calculate the total payload weight, using the number of passengers, crew and the weight of a person
    return (n_passengers * w_person + n_crew * w_person) * lbs_to_kg * g_0 + w_cargo * n_passengers * g_0


def calc_cruise_coefficient(cl, cd, r, velocity, c_j):
    # equations come from ADSEE-I lecture 2
    # print(r)
    # print(velocity)
    # print(cl/cd)
    return 1 / (np.e ** (r / ((velocity / (c_j * per_hr_to_N * g_0)) * (cl / cd))))


def calc_loiter_coefficient(aspect_ratio, oswald_factor, cd_0, e, c_j):
    CD = 2 * cd_0
    CL = np.sqrt((np.pi * aspect_ratio * oswald_factor * cd_0))

    return 1 / (np.e ** (e / ((1 / (c_j * g_0 * per_hr_to_N)) * (CL / CD))))


def calc_fuel_fraction(coefficients):
    # calculate the total required fuel by multiplying all phases of the flight profile
    fraction = 1
    for number in coefficients:
        fraction = fraction * number
    return 1 - fraction


def class_I(cl, cd, r_cruise, r_res, v_cruise, cj_cruise, W_tfo_frac, W_e_frac, fractions,
            N_pas, N_crew, W_person, W_cargo):  # , w_pay):
    cruise_1 = calc_cruise_coefficient(cl, cd, r_cruise, v_cruise, cj_cruise)
    # loiter = calc_loiter_coefficient(A, Oswald, S_ratio, C_fe, E, c_j_loiter)
    cruise_2 = calc_cruise_coefficient(cl, cd, r_res, v_cruise, cj_cruise)

    # print("The first cruise coefficient equals " + str(cruise_1))
    # print("The loiter coefficient equals " + str(loiter))
    # print("The second cruise coefficient equals " + str(cruise_2))

    mission_frac = np.array(
        [fractions[0], fractions[1], fractions[2], fractions[3], cruise_1, fractions[4], fractions[5], cruise_2,
         fractions[6], fractions[7]])
    # print("hoi")
    # print(W_e_frac)
    W_f_frac = calc_fuel_fraction(mission_frac)
    # print(W_f_frac)
    # print(W_e_frac)
    W_to_frac = 1 - W_f_frac - W_tfo_frac - W_e_frac

    W_P = calc_payload_weight(N_pas, N_crew, W_person, W_cargo)
    # W_P = w_pay
    W_TO = W_P / W_to_frac
    W_F = W_f_frac * W_TO
    W_E = W_e_frac * W_TO

    return np.array([W_TO, W_E, W_P, W_F])

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
