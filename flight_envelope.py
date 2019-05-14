import numpy as np
import matplotlib.pyplot as plt
from constants_and_conversions import *


def manoeuvring_envelope(w_to, h, cl_max_pos, s, v_cruise):
    # construct the manoeuvring plot
    n_max = min(3.8, max(2.5, 2.1 + (24000 / (w_to / lbs_to_kg + 10000))))
    rho = Rho_0 * ((1 + (a * h) / Temp_0) ** (-(g_0 / (R_gas * a) + 1)))

    cn_max_pos = 1.1 * cl_max_pos
    cn_max_neg = 0.8*cn_max_pos

    n_pos = np.arange(0., n_max + 0.1, 0.1)
    n_pos = np.append(n_pos, n_max)
    n_neg = np.arange(0., -1, -0.1)
    v_pos = np.zeros(len(n_pos))
    v_neg = np.zeros(len(n_neg))

    for i in range(len(n_pos)):
        v_pos[i] = np.sqrt((2 * w_to * n_pos[i]) / (rho * cn_max_pos * s))
    for i in range(len(n_neg)):
        if n_neg[i] == 0.:
            v_neg[i] = 0.
        else:
            v_neg[i] = -1 * np.sqrt((2 * w_to * -1 * n_neg[i]) / (rho * cn_max_neg * s))

    V_A = v_pos[-2]
    V_C = v_cruise
    V_D = 1.25 * V_C
    V_H = v_neg[-1]
    V_S = np.sqrt((2 * w_to) / (rho * cn_max_pos * s))

    v_pos[-1] = V_C
    v_pos = np.append(v_pos, [V_D, V_D])
    n_pos = np.append(n_pos, [n_pos[-1], 0])
    v_neg = np.append(v_neg, [V_C, V_D])
    n_neg = np.append(n_neg, [n_neg[-1], 0])

    speeds = np.array([V_S, 0, V_C, V_D])

    return v_pos, n_pos, v_neg, n_neg, speeds


def gust_envelope(w, h, cl_alpha, s, c, v_cruise, v):
    # construct gust loading plot
    # first determine density at altitude
    rho = Rho_0 * ((1 + (a * h) / Temp_0) ** (-(g_0 / (R_gas * a) + 1)))

    mu_g = (2 * (w / s)) / (rho * c * cl_alpha * g_0)
    # print(mu_g)
    K_g = 0.88 * mu_g / (5.3 + mu_g)
    # print("the kg value is " + str(K_g))

    # Kg = gust alleviation coefficient (as function of GH) with
    # c =  mean geometric chord (m)

    # Interpolate for the gust speeds at the desired altitude
    if h / ft_to_m < 20000.:
        U_B = 66 * ft_to_m
        U_C = 50 * ft_to_m
        U_D = 25 * ft_to_m
    else:
        U_B = (84.67 - 0.000933 * h / ft_to_m) * ft_to_m
        U_C = (66.67 - 0.000833 * h / ft_to_m) * ft_to_m
        U_D = (33.34 - 0.000417 * h / ft_to_m) * ft_to_m

    v_gusts = np.array([0, U_B, U_C, U_D])
    # print(v[0])
    #  Calculate V_B based on the stall speed and load increase
    V_B = v[0] * np.sqrt(1 + ((cl_alpha * K_g * U_B * v_cruise) / (w / s)))
    v[0] = 0
    v[1] = V_B

    n_pos = np.zeros(len(v))
    n_neg = np.zeros(len(v))
    # print("the gust speeds are " + str(v_gusts))
    # print("the speeds are " + str(v))
    #  Calculate the load factor based on the gust speed that accompanies the aircraft speed
    for i in range(len(v)):
        n_pos[i] = 1 + (0.5 * Rho_0 * cl_alpha * v_gusts[i] * v[i] * K_g) / (w / s)
        n_neg[i] = 1 - (0.5 * Rho_0 * cl_alpha * v_gusts[i] * v[i] * K_g) / (w / s)

    n_pos = np.append(n_pos, n_neg[-1])
    v_pos = np.append(v, v[-1])
    v_neg = v

    return v_pos, n_pos, v_neg, n_neg


def construct_envelope():
    # Note: used values are only estimation and are definitely not correct!
    v_pos, n_pos, v_neg, n_neg, gust_speeds = manoeuvring_envelope(1592281, 11000, 1.6, 265, 295)
    v_gust_pos, n_gust_pos, v_gust_neg, n_gust_neg = gust_envelope(1592281, 11000, 3.8, 295, 8, 265, gust_speeds)

    plt.plot(v_pos, n_pos)
    plt.plot(v_neg, n_neg)
    plt.plot(v_gust_pos, n_gust_pos)
    plt.plot(v_gust_neg, n_gust_neg)

    plt.show()


# construct_envelope()
