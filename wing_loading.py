from inputs import *
import matplotlib.pyplot as plt

#  Create the wing loading for which all inputs have to be evaluated
WS = np.arange(0., 8000., 50)
TW = np.arange(0., 0.6, 0.05)


#  Calculate take-off requirements, equation from ADSEE-I, lecture 3
def calc_take_off(cl, sigma, top, ws):
    tw = np.zeros((len(cl), len(ws)))
    label = []

    for i in range(len(tw)):
        label.append("Take-off with CL = " + str(round(cl[i], 2)))
        for j in range(len(tw[0])):
            tw[i][j] = (ws[j] / top) * (1 / cl[i]) * (1 / sigma)

    return tw, label


# Calculate the cruise requirement, equation from ADSEE-I, lecture 3
def calc_cruise(rho_cruise, cd_0, velocity, ws, aspect_ratio, oswald_factor):
    tw = np.zeros((len(aspect_ratio), len(ws)))
    label = []

    for i in range(len(tw)):
        label.append("Cruise with A = " + str(round(aspect_ratio[i], 2)))
        for j in range(len(tw[0])):
            # Check for ws being equal to 0
            if ws[j] == 0:
                ws[j] = 1

            tw[i][j] = ((Rho_0 / rho_cruise) ** 0.75) * ((cd_0 * 0.5 * rho_cruise * (velocity ** 2)) / ws[j] + ws[j] / (
                    np.pi * aspect_ratio[i] * oswald_factor * 0.5 * rho_cruise * (velocity ** 2)))

    return tw, label


# Calculate the climb rate requirement, equation from ADSEE-I, lecture 3
def calc_climb_rate(aspect_ratio, oswald_factor, cd_0, climb_rate, rho_to, ws):
    tw = np.zeros((len(aspect_ratio), len(ws)))
    cl = np.sqrt(3 * np.pi * cd_0 * oswald_factor) * np.sqrt(aspect_ratio)
    cd = 4 * cd_0
    label = []

    for i in range(len(tw)):
        label.append("c with A = " + str(round(aspect_ratio[i], 2)))
        for j in range(len(tw[0])):
            # Check for ws being equal to 0
            if ws[j] == 0:
                ws[j] = 1

            tw[i][j] = climb_rate / (np.sqrt(ws[j]) * np.sqrt(2 / (rho_to * cl[i]) + cd / cl[i]))

    return tw, label


# Calculate the climb gradient requirement, equation form ADSEE-I, lecture 3
def calc_climb_gradient(cd_0, delta_cd, n_engines, climb_gradient, oswald_factor, delta_oswald, aspect_ratio,
                        wing_loading):
    cd_0 += delta_cd
    oswald_factor += delta_oswald
    tw = np.zeros((len(aspect_ratio), len(wing_loading)))
    label = []

    for i in range(len(aspect_ratio)):
        label.append("c/V with A = " + str(round(aspect_ratio[i], 2)))
        for j in range(len(tw[0])):
            tw[i][j] = ((n_engines / (n_engines - 1)) * climb_gradient + 2 * np.sqrt(
                cd_0 / (np.pi * oswald_factor * aspect_ratio[i])))

    return tw, label


# Calculate the climb gradient requirement, equation form ADSEE-I, lecture 3
def calc_stall_loading(v_stall, cl, rho, tw):
    ws = np.zeros((len(cl), len(tw)))

    for i in range(len(ws)):

        for j in range(len(ws[0])):
            ws[i][j] = 0.5 * rho[i] * cl[i] * (v_stall[i] ** 2)

    return ws, ["Stall in take-off config", "Stall in cruise config", "Stall in landing config"]


# Calculate the climb gradient requirement, equation form ADSEE-I, lecture 3
def calc_landing_loading(f, rho, v_stall, cl, tw):
    ws = np.zeros((len(cl), len(tw)))
    label = []

    for i in range(len(ws)):
        label.append("Landing with CL_Land_max = " + str(round(cl[i], 2)))
        for j in range(len(ws[0])):
            ws[i][j] = (0.5 * rho * cl[i] * (v_stall ** 2)) / f

    return ws, label


thrust_to_weight_to, labels_to = calc_take_off(CL_TO, Sigma_TO, TOP, WS)
thrust_to_weight_cruise, labels_cruise = calc_cruise(Rho_Cruise, CD_0, V_cruise, WS, AR, Oswald_factor)
thrust_to_weight_climb_rate, labels_climb_rate = calc_climb_rate(AR, Oswald_factor, CD_0, Climb_rate, Rho_TO, WS)
thrust_to_weight_climb_gradient, labels_climb_gradient = calc_climb_gradient(CD_0, Delta_CD_TO_gear_up, N_engines,
                                                                             Climb_gradient,
                                                                             Oswald_factor, Delta_oswald_TO, AR, WS)

wing_loading_landing, labels_landing = calc_landing_loading(Landing_factor, Rho_Landing, V_stall_Landing, CL_Landing,
                                                            TW)
wing_loading_stall, labels_stall = calc_stall_loading([V_stall_TO, V_stall_Cruise, V_stall_Landing],
                                                      [CL_TO_max, CL_Cruise_max, CL_Landing_max],
                                                      [Rho_TO, Rho_Cruise, Rho_Landing],
                                                      TW)

# print("The take-off requirement is seb by :" + str(thrust_to_weight_to))
# print("The cruise requirement is set by :" + str(thrust_to_weight_cruise))
# print("The climb rate requirement is set by :" + str(thrust_to_weight_climb_rate))
# print("the climb gradient requirement is set by :" + str(thrust_to_weight_climb_gradient))
# print("the stall speed requirement is set by :" + str(wing_loading_stall))
# print("the landing requirement is set by :" + str(wing_loading_landing))

all_requirements = [thrust_to_weight_climb_gradient, thrust_to_weight_climb_rate, thrust_to_weight_cruise,
                    thrust_to_weight_to, wing_loading_landing, wing_loading_stall]
labels = [labels_climb_gradient, labels_climb_rate, labels_cruise, labels_to, labels_landing, labels_stall]

for p in range(len(all_requirements)):
    for q in range(len(all_requirements[p])):
        if len(all_requirements[p][q]) == len(TW):
            plt.plot(all_requirements[p][q], TW, label=labels[p][q])
        else:
            plt.plot(WS, all_requirements[p][q], label=labels[p][q])

plt.legend(loc='upper right')
plt.ylim(0.0, 0.4)
plt.show()
