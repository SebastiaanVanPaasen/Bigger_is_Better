from inputs import *


def calc_payload_weight(n_passengers, n_crew, w_person):
    # calculate the total payload weight, using the number of passengers, crew and the weight of a person
    return (n_passengers * w_person * lbs_to_kg + n_crew * w_person * lbs_to_kg) * g_0


def calc_cruise_coefficient(aspect_ratio, oswald_factor, surface_ratio, c_fe, r, velocity, c_j):
    # c_fe is the equivalent skin friction coeficient
    # surface_ratio is the wetted surface over the total surface
    # equations come from ADSEE-I lecture 2

    cd_0 = c_fe * surface_ratio
    LD_ratio = 0.75 * np.sqrt((np.pi * aspect_ratio * oswald_factor) / (3 * cd_0))

    return 1 / (np.e ** (r / ((velocity / (c_j * per_hr_to_N * g_0)) * LD_ratio)))


def calc_loiter_coefficient(aspect_ratio, oswald_factor, surface_ratio, c_fe, e, c_j):
    CD_0 = c_fe * surface_ratio
    CD = 2 * CD_0
    CL = np.sqrt((np.pi * aspect_ratio * oswald_factor * CD_0))

    return 1 / (np.e ** (e / ((1 / (c_j * g_0 * per_hr_to_N)) * (CL / CD))))


def calc_fuel_fraction(coefficients):
    # calculate the total required fuel by multiplying all phases of the flight profile
    fraction = 1
    for number in coefficients:
        fraction = fraction * number
    return 1 - fraction


def class_I(aspect_ratio, oswald, s_ratio, c_fe, r_cruise, r_res, v_cruise, cj_cruise):
    cruise_1 = calc_cruise_coefficient(aspect_ratio, oswald, s_ratio, c_fe, r_cruise, v_cruise, cj_cruise)
    # loiter = calc_loiter_coefficient(A, Oswald, S_ratio, C_fe, E, c_j_loiter)
    cruise_2 = calc_cruise_coefficient(aspect_ratio, oswald, s_ratio, c_fe, r_res, v_cruise, cj_cruise)

    # print("The first cruise coefficient equals " + str(cruise_1))
    # print("The loiter coefficient equals " + str(loiter))
    # print("The second cruise coefficient equals " + str(cruise_2))

    mission_frac = np.array([start, taxi, t_o, climb_1, cruise_1, descent_1, climb_2, cruise_2, descent_2, landing])

    W_f_frac = calc_fuel_fraction(mission_frac)
    W_to_frac = 1 - W_f_frac - W_tfo_frac - W_e_frac

    # print(w_f_frac)
    # print(w_to_frac)

    W_P = calc_payload_weight(N_pas, N_crew, W_person)
    W_TO = W_P / W_to_frac
    W_F = W_f_frac * W_TO
    W_E = W_e_frac * W_TO

    return np.array([W_TO, W_E, W_P, W_F]) / g_0
