import numpy as np
import matplotlib.pyplot as plt

lb_to_kg = 0.45359237
hr_to_sec = 3600
kts_to_ms = 0.51444444444
nm_to_km = 1.852
g_0 = 9.80565  # m/s^2
per_hr_to_N = 1 / (hr_to_sec * g_0)


def calc_payload_weight(n_passengers, n_crew, w_person):
    return (n_passengers * w_person * lb_to_kg + n_crew * w_person * lb_to_kg) * g_0


def calc_cruise_coefficient(aspect_ratio, oswald_factor, surface_ratio, c_fe, r, velocity, c_j):
    cd_0 = c_fe * surface_ratio
    LD_ratio = 0.75 * np.sqrt((np.pi * aspect_ratio * oswald_factor) / (3 * cd_0))

    return 1 / (np.e ** (r / ((velocity / (c_j * per_hr_to_N * g_0)) * LD_ratio)))


# def calc_loiter_coefficient(aspect_ratio, oswald_factor, surface_ratio, c_fe, e, c_j):
#     CD_0 = c_fe * surface_ratio
#     CD = 2 * CD_0
#     CL = np.sqrt((np.pi * aspect_ratio * oswald_factor * CD_0))
#
#     return 1 / (np.e ** (e / ((1 / (c_j * g_0 * per_hr_to_N)) * (CL / CD))))


def calc_fuel_fraction(coefficients):
    fraction = 1
    for number in coefficients:
        fraction = fraction * number
    return 1 - fraction


# ----------------------------------------------------------------------------------------------------------------------
# Define fuel fractions from statistics
start = 0.99
taxi = 0.99
t_o = 0.995
climb_1 = 0.98
descent_1 = 0.99
climb_2 = 0.98
descent_2 = 0.99
landing = 0.992
# ----------------------------------------------------------------------------------------------------------------------
N_pas = 450
N_crew = 11
W_person = 205  # lbs
W_tfo = 0.003  # estimated from slides ADSEE-I, lecture 3
empty_fraction = 0.525  # Based on Ed Obart
# max_fuel_fraction = 0.410283406
# ----------------------------------------------------------------------------------------------------------------------


def main(r_1, r_2):
    V = 499 * kts_to_ms  # m/s
    C_fe = 0.003  # estimated
    S_ratio = 6.3  # estimated
    A = 6.96  # estimated
    oswald = 0.8  # estimated
    c_j_cruise = 0.75  # 1/hr
    c_j_loiter = 0.5  # 1/hr

    cruise_1 = calc_cruise_coefficient(A, oswald, S_ratio, C_fe, r_1, V, c_j_cruise)
    # loiter = calc_loiter_coefficient(A, oswald, S_ratio, C_fe, E, c_j_loiter)
    cruise_2 = calc_cruise_coefficient(A, oswald, S_ratio, C_fe, r_2, V, c_j_cruise)

    print("The first cruise coefficient equals " + str(cruise_1))
    # print("The loiter coefficient equals " + str(loiter))
    print("The second cruise coefficient equals " + str(cruise_2))

    fuel_fractions = np.array(
        [start, taxi, t_o, climb_1, cruise_1, descent_1, climb_2, cruise_2, descent_2, landing])

    fuel_fraction = calc_fuel_fraction(fuel_fractions)
    print(fuel_fraction)
    weight_fraction = 1 - fuel_fraction - W_tfo - empty_fraction
    print(weight_fraction)

    W_P = calc_payload_weight(N_pas, N_crew, W_person)
    W_TO = W_P / weight_fraction
    W_F = fuel_fraction * W_TO
    W_E = empty_fraction * W_TO

    return np.array([W_TO, W_E, W_P, W_F]) / g_0


# R_1 = np.arange(1000000., 2100000., 100000.)
R_2 = 250 * 1.852 * 1000  # m
# results = np.zeros((4, len(R_1)))
#
# for i in range(len(R_1)):
#     value = main(R_1[i], R_2)
#     for j in range(4):
#         results[j][i] = value[j]

# print((design_range[1] - different_range[1]) / ((R_1 - D_range) / 1000))
design_range = main(1800000, R_2)
print("The take-off mass equals " + str(round(design_range[0], 2)))
print("The fuel mass equals " + str(round(design_range[3], 2)))
print("The payload mass equals " + str(round(design_range[2], 2)))
print("The empty mass equals " + str(round(design_range[1], 2)))
#
# print("The take-off mass equals " + str(round(different_range[1], 2)))
# print("The fuel mass equals " + str(round(different_range[2], 2)))
# print("The payload mass equals " + str(round(different_range[0], 2)))
# print("The empty mass equals " + str(round(different_range[3], 2)))
