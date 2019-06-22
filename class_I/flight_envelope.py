import numpy as np
import matplotlib.pyplot as plt
import constants_and_conversions as cc


def manoeuvring_envelope(w_to, h, cl_max_pos, S, v_cruise):
    # construct the manoeuvring plot
    n_max = min(3.8, max(2.5, 2.1 + (24000 / (((w_to / cc.g_0) / cc.lbs_to_kg) + 10000))))
    n_min = -1
    
    rho = cc.Rho_0 * ((1 + (cc.a * h) / cc.Temp_0) ** (-(cc.g_0 / (cc.R_gas * cc.a) + 1)))
    cn_max_pos = 1.1*cl_max_pos
    print(rho)
    cn_max_neg =1.1* 0.8*cl_max_pos#0.8 * cn_max_pos

    n_lim_pos = np.arange(0., n_max + 0.1, 0.1)
    n_lim_pos = np.append(n_lim_pos, n_max)
    n_lim_neg = np.arange(0., -1 - 0.1, -0.1)
#    print(n_lim_neg)
    
    v_pos = np.zeros(len(n_lim_pos))
    v_neg = np.zeros(len(n_lim_neg))

    for i in range(len(n_lim_pos)):
        if n_lim_pos[i] == 0.:
            v_pos[i] = 0.
        else:
            v_pos[i] = np.sqrt((2 * w_to * n_lim_pos[i]) / (rho * cn_max_pos * S))
        # print(v_pos[i])
        
    for i in range(len(n_lim_neg)):
        if n_lim_neg[i] == 0.:
            v_neg[i] = 0.
        else:
            v_neg[i] = np.sqrt((2 * w_to * -1 * n_lim_neg[i]) / (rho * cn_max_neg * S))
    
    V_C = v_cruise
    V_D = 1.25 * V_C
    V_S = np.sqrt((2 * w_to) / (rho * cn_max_pos * S))
    V_A = V_S * np.sqrt(n_max)

    v_pos = np.append(v_pos, [V_D, V_D])
    n_lim_pos = np.append(n_lim_pos, [n_lim_pos[-1], 0])
    v_neg = np.append(v_neg, [V_C, V_D])
    n_lim_neg = np.append(n_lim_neg, [-1, 0])

    speeds = np.array([0, V_S, V_A, V_C, V_D])
    
    
#    print("hoi", v_pos[-2], V_A)
#    print(speeds)
    
    
    return v_pos,v_neg, n_lim_pos, n_lim_neg, speeds

#Vs1 = manoeuvring_envelope(1828542.22, 10000, 1.6 , 0.1 ,280, 239.28)[4][1] 
#print(Vs1)
def gust_envelope(w_to, h, cl_alpha, S, c, v_cruise, Vs1):
    # construct gust loading plot
    # first determine density at altitude
    rho = cc.Rho_0 * ((1 + (cc.a * h) / cc.Temp_0) ** (-(cc.g_0 / (cc.R_gas * cc.a) + 1)))
    rho = rho*cc.kgm3_to_slugft3
    
    mu_g = (2 * (w_to / S)/cc.psf_to_nm2) / (rho * (c / cc.ft_to_m) * cl_alpha * (3.2808399 * cc.g_0))
#    print(mu_g)
    K_g = 0.88 * mu_g / (5.3 + mu_g)
#    print("the kg value is " + str(K_g))

    # Kg = gust alleviation coefficient (as function of GH) with
    # c =  mean geometric chord (m)

    # Interpolate for the gust speeds at the desired altitude
    if h / cc.ft_to_m < 20000.:
        U_B = 66 * cc.ft_to_m
        U_C = 50 * cc.ft_to_m
        U_D = 25 * cc.ft_to_m
    else:
        U_B = (84.67 - 0.000933 * h / cc.ft_to_m) # * cc.ft_to_m
        U_C = (66.67 - 0.000833 * h / cc.ft_to_m) #* cc.ft_to_m
        U_D = (33.34 - 0.000417 * h / cc.ft_to_m) #* cc.ft_to_m

    
    v_gusts = np.array([0, U_B, U_C, U_D])
    # print(v[0])
    #  Calculate V_B based on the stall speed and load increase
    # V_B = v[0] * np.sqrt(1 + ((cl_alpha * K_g * U_B * v_cruise) / (w / s)))
    V_C = v_cruise
    V_D = 1.25 * V_C
    V_B_higher = (Vs1 / cc.kts_to_ms) * np.sqrt(1 + K_g * 44 * (V_C / cc.kts_to_ms) * cl_alpha/(498*(w_to / S)/cc.psf_to_nm2))
    V_B_higher = V_B_higher* cc.kts_to_ms
    V_B_lower = v_cruise - 43 * cc.kts_to_ms # still needs to be changed 
    V_B = (V_B_higher + V_B_lower)/2
#    print(V_C)
#    print("vb",V_B_higher)
#    print("vblow",V_B_lower)
    V = [0, V_B, V_C, V_D]
    n_lim_pos = np.zeros(len(V))
    n_lim_neg = np.zeros(len(V))
    # print("the gust speeds are " + str(v_gusts))
    # print("the speeds are " + str(v))
    #  Calculate the load factor based on the gust speed that accompanies the aircraft speed
    for i in range(len(V)):
        n_lim_pos[i] = 1 + (cl_alpha * v_gusts[i] * V[i]/cc.kts_to_ms * K_g) / (498*(w_to / S)/cc.psf_to_nm2)
        n_lim_neg[i] = 1 - (cl_alpha * v_gusts[i] * V[i]/cc.kts_to_ms * K_g) / (498*(w_to / S)/cc.psf_to_nm2)
#    n_pos = np.append(n_pos, n_neg[-1])
#    v_pos = np.append(v, v[-1])
#    v_neg = v
#    print(n_lim_pos)
    return V, n_lim_pos, n_lim_neg


def construct_envelope():
    # Note: used values are only estimation and are definitely not correct!
    v_pos,v_neg, n_lim_pos, n_lim_neg,speeds = manoeuvring_envelope(1520276, 9000, 1.52, 211.2, 218)
    V_gust, n_gust_pos, n_gust_neg = gust_envelope(1520276, 9000, 6.2, 211, 4.03, 218, speeds[1]) #w_to, h, cl_alpha, S, c, v_cruise, Vs1
#    print(n_gust_pos, n_gust_neg)
    plt.plot(v_pos, n_lim_pos)
    plt.plot(v_neg, n_lim_neg)
    plt.plot(V_gust,n_gust_pos)
    plt.plot(V_gust,n_gust_neg)
    
    plt.xlabel("V [m/s]")
    plt.ylabel("n [-]")

    plt.show()

    return V_gust, n_gust_pos, n_gust_neg, speeds


print(construct_envelope())

#print(construct_envelope())
