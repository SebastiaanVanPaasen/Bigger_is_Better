import constants_and_conversions as cc
import numpy as np
import matplotlib.pyplot as plt


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
            label.append("Take-off with CL_To_max = " + str(round(cl[i], 2)))
            for j in range(len(tw[0])):
                tw[i][j] = (self.ws[j] / top) * (1 / cl[i]) * (1 / self.sigma)

        return tw, label

    # Calculate the cruise requirement, equation from ADSEE-I, lecture 3
    def calc_cruise(self, rho_cruise, velocity):
        tw = np.zeros((len(self.aspect_ratio), len(self.ws)))
        label = []

        for i in range(len(tw)):
            label.append("Cruise config")  # = " + str(round(self.aspect_ratio[i], 2)))
            for j in range(len(tw[0])):
                # Check for ws being equal to 0
                if self.ws[j] == 0:
                    self.ws[j] = 1

                tw[i][j] = ((cc.Rho_0 / rho_cruise) ** 0.75) * (
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
            label.append("Climb rate ")  # + str(round(self.aspect_ratio[i], 2)))
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
            label.append("Climb gradient ")  # + str(round(self.aspect_ratio[i], 2)))
            for j in range(len(tw[0])):
                tw[i][j] = ((n_engines / (n_engines - 1)) * (climb_gradient + 2 * np.sqrt(
                    cd_0 / (np.pi * oswald_factor * self.aspect_ratio[i]))))

        return tw, label

    # Calculate the climb gradient requirement, equation form ADSEE-I, lecture 3
    def calc_stall(self, v_stall, cl, rho):
        ws = np.zeros((len(cl), len(self.tw)))

        for i in range(len(ws)):
            for j in range(len(ws[0])):
                ws[i][j] = 0.5 * rho[i] * cl[i] * (v_stall[i] ** 2)

        return ws, ["Stall in cruise config", "Stall in landing config"]

    # Calculate the climb gradient requirement, equation form ADSEE-I, lecture 3
    def calc_landing(self, f, rho, v_stall, cl):
        ws = np.zeros((len(cl), len(self.tw)))
        label = []

        for i in range(len(ws)):
            label.append("Landing with CL_Land_max = " + str(round(cl[i], 2)))
            for j in range(len(ws[0])):
                ws[i][j] = (0.5 * rho * cl[i] * (v_stall ** 2)) / f

        return ws, label


def plot_diagram(all_requirements, labels, ws, tw):
    r = 0
    plt.figure(figsize=(15, 9))
    for p in range(len(all_requirements)):
        for q in range(len(all_requirements[p])):
            if len(all_requirements[p][q]) == len(tw):
                plt.plot(all_requirements[p][q], tw, label=labels[p][q], marker="o")
            else:
                if r < 10:
                    plt.plot(ws, all_requirements[p][q], label=labels[p][q], marker="^")
                else:
                    plt.plot(ws, all_requirements[p][q], label=labels[p][q], marker="o")
                r += 1

    plt.scatter(5860, 0.2819, label="Airbus A330-300 (wb)", c="r", marker="o", s=100)
    plt.scatter(6940, 0.2205, label="Airbus A340-200 (wb)", c="g", marker="o", s=100)
    plt.scatter(7318, 0.2272, label="Airbus A340-300 (wb)", c="b", marker="o", s=100)
    plt.scatter(8184, 0.2634, label="Airbus A340-500 (wb)", c="y", marker="o", s=100)
    plt.scatter(8184, 0.2783, label="Airbus A340-600 (wb)", c="k", marker="o", s=100)
    plt.scatter(5562, 0.2878, label="Boeing 777-200 (wb)", c="m", marker="o", s=100)
    plt.scatter(6576, 0.2655, label="Boeing 777-200IGW (wb)", c="c", marker="o", s=100)
    plt.scatter(6861, 0.2818, label="Boeing 777-200XI (wb)", c="r", marker="^", s=100)
    plt.scatter(7173, 0.2841, label="Boeing 777-200X2 (wb)", c="g", marker="^", s=100)
    plt.scatter(5479, 0.322, label="Tupolev Tu-334 (nb)", c="b", marker="^", s=100)
    plt.scatter(5178, 0.272, label="BAe RJ70 (nb)", c="y", marker="^", s=100)
    plt.scatter(5351, 0.301, label="BAe RJ85 (nb)", c="k", marker="^", s=100)
    plt.scatter(5610, 0.287, label="BAe RJ100 (nb)", c="m", marker="^", s=100)
    plt.scatter(5840, 0.276, label="BAe RJ115 (nb)", c="c", marker="^", s=100)
    plt.scatter(3678, 0.333, label="Embraer EMB-145 (nb)", c="r", marker="*", s=100)
    plt.scatter(3853, 0.342, label="Fokker F70 (nb)", c="g", marker="*", s=100)
    plt.scatter(4519, 0.291, label="Fokker F100 (nb)", c="b", marker="*", s=100)
    plt.scatter(6000, 0.23, label="Design point minimum T/W", c="r", marker="v", s=100)
    plt.scatter(7000, 0.26, label="Design point average T/W and W/S", c="y", marker="v", s=100)
    plt.scatter(8000, 0.29, label="Design point maximum W/S", c="b", marker="v", s=100)

    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel(r"$\frac{W}{S}$ [$N/m^2$]")
    plt.ylabel(r"$\frac{T}{W}$ [-]")
    plt.xlim(0, 10000)
    plt.ylim(0.0, 0.4)
    plt.show()
#    plt.savefig("Wing_loading_diagram")

# Landing requirements -------------------------------------------------------------------------------------------------
Landing_runway = 2500.  # m  Guestimated from ADSEE-I L3
Landing_factor = 0.84  # Guestimated from ADSEE-I L3

# Take-off paramters ---------------------------------------------------------------------------------------------------
TOP = 220. * cc.lbft2_Nm2  # Guestimated from ADSEE-I L3
Sigma_TO = 1  # Ratio of densities

# Stall speeds ---------------------------------------------------------------------------------------------------------
V_stall_Cruise = 228 / 1.2  # m/s   Guestimated from ADSEE-I L3 Take requirements
V_stall_Landing = min(np.sqrt(Landing_runway / 0.5847), 65.)  # m/s     Guestimated from ADSEE-I L3

# Lift and drag coefficients -------------------------------------------------------------------------------------------
CL_Cruise_max = 1.3  # Guestimated from ADSEE-I L3
CL_Landing_max = 3.2  # From Obert
CL_TO_max = 0.9 * CL_Landing_max  # Guestimated from ADSEE-I L3

# Changes in coefficients due to landing or take-off -------------------------------------------------------------------
Delta_CD0_TO_gear_up = 0.015  # Guestimated from ADSEE-I L3
Delta_CD0_TO_gear_down = Delta_CD0_TO_gear_up + 0.02  # Guestimated from ADSEE-I L3
Delta_CD0_Land = 0.085  # Guestimated from ADSEE-I L3

Delta_oswald_TO = 0.05  # Guestimated from ADSEE-I L3
Delta_oswald_Land = 0.1  # Guestimated from ADSEE-I L3

# Estimated parameters -------------------------------------------------------------------------------------------------
Climb_rate = 2000. * cc.ftmin_to_ms  # m/s      set climb rate by CS25
Climb_gradient = 0.024  # c/V set by CS25 in case 4 engines, should be 0.03

# Create ranges
CL_TO = [CL_TO_max]  # / 1.21 - 0.2), CL_TO_max / 1.21, (CL_TO_max / 1.21 + 0.2)]  # Based on CL_max
CL_Landing = [CL_Landing_max]  # - 0.2), CL_Landing_max, (CL_Landing_max + 0.2)]  # Based on CL_max
AR = [14]  # - 2, A, A + 2]  # Based on A in class I
WS = np.arange(0., 10250., 250)
TW = np.arange(0., 0.45, 0.05)

H_cr = 9000
V_cr = 228
Rho_cr = cc.Rho_0 * ((1 + (cc.a * H_cr) / cc.Temp_0) ** (-(cc.g_0 / (cc.R_gas * cc.a))))
N_engines = 2

def final_diagram(cd0, os):
    diagram = TWdiagram(Sigma_TO, WS, TW, cd0, os, AR)
    take_off, labels_to = diagram.calc_take_off(CL_TO, TOP)
    cruise, labels_cruise = diagram.calc_cruise(Rho_cr, 220)
    climb_rate_performance, labels_climb_rate = diagram.calc_climb_rate(Climb_rate, cc.Rho_0)
    climb_gradient_performance, labels_climb_gradient = diagram.calc_climb_gradient(Delta_CD0_TO_gear_up, N_engines,
                                                                                    Climb_gradient, Delta_oswald_TO)
    landing, labels_landing = diagram.calc_landing(Landing_factor, cc.Rho_0, V_stall_Landing, CL_Landing, )
    stall, labels_stall = diagram.calc_stall([V_stall_Cruise, V_stall_Landing], [CL_Cruise_max, CL_Landing_max],
                                             [Rho_cr, cc.Rho_0])

    requirement = [climb_gradient_performance, climb_rate_performance, cruise, take_off, landing, stall]
    labelling = [labels_climb_gradient, labels_climb_rate, labels_cruise, labels_to, labels_landing, labels_stall]

    plot_diagram(requirement, labelling, WS, TW)


final_diagram(0.02589, 0.85)
