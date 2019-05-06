# Put all the inputs that you use in other files in here, the file has been divided in different sections for constants
# input parameters and conversion factors
import numpy as np

# Standard constants ---------------------------------------------------------------------------------------------------
Rho_0 = 1.225  # kg/m^3
a = -0.065  # K/m up to 10 km

# Conversion factors ---------------------------------------------------------------------------------------------------
ftmin_to_ms = 0.00508
kts_to_ms = 0.51444444

# Input parameters for the aircraft ------------------------------------------------------------------------------------
Landing_runway = 2500  # m  guestimated
Landing_factor = 0.84  # guestimated

Rho_TO = 1.225  # kg/m^3    standard sea-level density
Rho_Cruise = 0.6  # kg/m^3      estimated cruise density
Rho_Landing = 1.225  # kg/m^3   standard sea-level density

TOP = 4325  # guestimated
Sigma_TO = Rho_TO / Rho_0  # ratio of densities

V_stall_TO = 60  # m/s  guestimated
V_stall_Cruise = 180 * kts_to_ms  # m/s   guestimated
V_stall_Landing = np.sqrt(Landing_runway / 0.5847)  # m/s     based on runway length

CL_TO_max = 2.2  # guestimated
CL_Cruise_max = 1.6  # guestimated
CL_Landing_max = 2.55  # guestimated

CL_TO = [(CL_TO_max / 1.21 - 0.2), CL_TO_max / 1.21, (CL_TO_max / 1.21 + 0.2)]  # Estimated possible cl values
CL_Landing = [(CL_Landing_max - 0.2), CL_Landing_max, (CL_Landing_max + 0.2)]  # guestimated


CD_0 = 0.018  # guestimated
AR = [4.96, 6.96, 8.96]  # Estimated possible aspect ratio values
V_cruise = 250.  # m/s      estimated cruise speed
Oswald_factor = 0.76  # guestimated

Climb_rate = 2000 * ftmin_to_ms  # m/s      set climb rate by CS25
Climb_gradient = 0.03  # c/V set by CS25
N_engines = 4  # guestimated

Delta_CD_TO_gear_up = 0.015  # estimated from slides ADSEE-I, lecture 3
Delta_CD_TO_gear_down = Delta_CD_TO_gear_up + 0.02  # estimated from slides ADSEE-I, lecture 3
Delta_CD_Land = 0.085  # estimated from slides ADSEE-I, lecture 3

Delta_oswald_TO = 0.05  # estimated from slides ADSEE-I, lecture 3
Delta_oswald_Land = 0.1  # estimated from slides ADSEE-I, lecture 3
