from inputs import *
import matplotlib.pyplot as plt

#  Create the wing loading for which all inputs have to be evaluated
WS = np.arange(0., 8050., 50)
TW = np.arange(0., 0.45, 0.05)


class TWdiagram:
    def __init__(self, sigma, ws, tw, cd_0, oswald_factor, aspect_ratio):
        self.sigma = sigma
        self.ws = ws
        self.tw = tw
        self.cd_0 = cd_0
        self.oswald = oswald_factor
        self.aspect_ratio = aspect_ratio

    #  Calculate take-off requirements, equation from ADSEE-I, lecture 3
    def calc_take_off(self, cl, top):
        tw = np.zeros((len(cl), len(self.ws)))
        label = []

        for i in range(len(tw)):
            label.append("Take-off with CL = " + str(round(cl[i], 2)))
            for j in range(len(tw[0])):
                tw[i][j] = (self.ws[j] / top) * (1 / cl[i]) * (1 / self.sigma)

        return tw, label

    # Calculate the cruise requirement, equation from ADSEE-I, lecture 3
    def calc_cruise(self, rho_cruise, velocity):
        tw = np.zeros((len(self.aspect_ratio), len(self.ws)))
        label = []

        for i in range(len(tw)):
            label.append("Cruise with A = " + str(round(self.aspect_ratio[i], 2)))
            for j in range(len(tw[0])):
                # Check for ws being equal to 0
                if self.ws[j] == 0:
                    self.ws[j] = 1

                tw[i][j] = ((Rho_0 / rho_cruise) ** 0.75) * (
                        (self.cd_0 * 0.5 * rho_cruise * (velocity ** 2)) / self.ws[j] + self.ws[j] / (
                            np.pi * self.aspect_ratio[i] * self.oswald * 0.5 * rho_cruise * (velocity ** 2)))

        return tw, label

    # Calculate the climb rate requirement, equation from ADSEE-I, lecture 3
    def calc_climb_rate(self, climb_rate, rho_to):
        tw = np.zeros((len(self.aspect_ratio), len(self.ws)))
        cl = np.sqrt(3 * np.pi * self.cd_0 * self.oswald) * np.sqrt(self.aspect_ratio)
        cd = 4 * self.cd_0
        label = []

        for i in range(len(tw)):
            label.append("c with A = " + str(round(self.aspect_ratio[i], 2)))
            for j in range(len(tw[0])):
                # Check for ws being equal to 0
                if self.ws[j] == 0:
                    self.ws[j] = 1

                tw[i][j] = climb_rate / (np.sqrt(self.ws[j]) * np.sqrt(2 / (rho_to * cl[i]) + cd / cl[i]))

        return tw, label

    # Calculate the climb gradient requirement, equation form ADSEE-I, lecture 3
    def calc_climb_gradient(self, delta_cd, n_engines, climb_gradient, delta_oswald):
        cd_0 = self.cd_0 + delta_cd
        oswald_factor = self.oswald + delta_oswald
        tw = np.zeros((len(self.aspect_ratio), len(self.ws)))
        label = []

        for i in range(len(self.aspect_ratio)):
            label.append("c/V with A = " + str(round(self.aspect_ratio[i], 2)))
            for j in range(len(tw[0])):
                tw[i][j] = ((n_engines / (n_engines - 1)) * climb_gradient + 2 * np.sqrt(
                    cd_0 / (np.pi * oswald_factor * self.aspect_ratio[i])))

        return tw, label

    # Calculate the climb gradient requirement, equation form ADSEE-I, lecture 3
    def calc_stall(self, v_stall, cl, rho):
        ws = np.zeros((len(cl), len(self.tw)))

        for i in range(len(ws)):

            for j in range(len(ws[0])):
                ws[i][j] = 0.5 * rho[i] * cl[i] * (v_stall[i] ** 2)

        return ws, ["Stall in take-off config", "Stall in cruise config", "Stall in landing config"]

    # Calculate the climb gradient requirement, equation form ADSEE-I, lecture 3
    def calc_landing(self, f, rho, v_stall, cl):
        ws = np.zeros((len(cl), len(self.tw)))
        label = []

        for i in range(len(ws)):
            label.append("Landing with CL_Land_max = " + str(round(cl[i], 2)))
            for j in range(len(ws[0])):
                ws[i][j] = (0.5 * rho * cl[i] * (v_stall ** 2)) / f

        return ws, label


diagram = TWdiagram(Sigma_TO, WS, TW, CD_0, Oswald_factor, AR)
take_off, labels_to = diagram.calc_take_off(CL_TO, TOP)
cruise, labels_cruise = diagram.calc_cruise(Rho_Cruise, V_cruise)
climb_rate_performance, labels_climb_rate = diagram.calc_climb_rate(Climb_rate, Rho_TO)
climb_gradient_performance, labels_climb_gradient = diagram.calc_climb_gradient(Delta_CD_TO_gear_up, N_engines,
                                                                                Climb_gradient, Delta_oswald_TO)

landing, labels_landing = diagram.calc_landing(Landing_factor, Rho_Landing, V_stall_Landing, CL_Landing, )
stall, labels_stall = diagram.calc_stall([V_stall_TO, V_stall_Cruise, V_stall_Landing],
                                         [CL_TO_max, CL_Cruise_max, CL_Landing_max],
                                         [Rho_TO, Rho_Cruise, Rho_Landing],
                                         )

# print("The take-off requirement is seb by :" + str(thrust_to_weight_to))
# print("The cruise requirement is set by :" + str(thrust_to_weight_cruise))
# print("The climb rate requirement is set by :" + str(thrust_to_weight_climb_rate))
# print("the climb gradient requirement is set by :" + str(thrust_to_weight_climb_gradient))
# print("the stall speed requirement is set by :" + str(wing_loading_stall))
# print("the landing requirement is set by :" + str(wing_loading_landing))

all_requirements = [climb_gradient_performance, climb_rate_performance, cruise,
                    take_off, landing, stall]
labels = [labels_climb_gradient, labels_climb_rate, labels_cruise, labels_to, labels_landing, labels_stall]

for p in range(len(all_requirements)):
    for q in range(len(all_requirements[p])):
        if len(all_requirements[p][q]) == len(TW):
            plt.plot(all_requirements[p][q], TW, label=labels[p][q])
        else:
            plt.plot(WS, all_requirements[p][q], label=labels[p][q])

plt.legend(loc='upper right')
plt.xlabel("W/S [N/m^2]")
plt.ylabel("T/W [-]")
plt.xlim(0, 8000)
plt.ylim(0.0, 0.4)
plt.show()
