# Put all the inputs that you use in other files in here, the file has been divided in different sections for constants
# input parameters and conversion factors
import numpy as np

# Standard constants ---------------------------------------------------------------------------------------------------
Rho_0 = 1.225  # kg/m^3
a = -0.065  # K/m up to 10 km

# Conversion factors ---------------------------------------------------------------------------------------------------
ftmin_to_ms = 0.00508
kts_to_ms = 0.51444444
lbft2_Nm2 = 47.880172

# Input parameters for the aircraft ------------------------------------------------------------------------------------
Landing_runway = 2500  # m  Guestimated from ADSEE-I L3
Landing_factor = 0.84  # Guestimated from ADSEE-I L3

Rho_TO = 1.225  # kg/m^3    standard sea-level density
Rho_Cruise = 0.6  # kg/m^3      estimated cruise density
Rho_Landing = 1.225  # kg/m^3   standard sea-level density

TOP = 220 * lbft2_Nm2  # Guestimated from ADSEE-I L3
Sigma_TO = Rho_TO / Rho_0  # Ratio of densities

V_stall_Cruise = 140  # m/s   Guestimated from ADSEE-I L3
V_stall_Landing = min(np.sqrt(Landing_runway / 0.5847), 65)  # m/s     Guestimated from ADSEE-I L3
V_stall_TO = V_stall_Landing * 1.05  # m/s  guestimated

CL_Cruise_max = 1.6  # Guestimated from ADSEE-I L3
CL_Landing_max = 3.2  # Guestimated from ADSEE-I L3
CL_TO_max = 0.8 * CL_Landing_max  # Guestimated from ADSEE-I L3

CL_TO = [(CL_TO_max / 1.21 - 0.2), CL_TO_max / 1.21, (CL_TO_max / 1.21 + 0.2)]  # Based on CL_max
CL_Landing = [(CL_Landing_max - 0.2), CL_Landing_max, (CL_Landing_max + 0.2)]  # Based on CL_max

CD_0 = 0.018  # guestimated from ADSEE-I L3
AR = [4.96, 6.96, 8.96]  #
V_cruise = 250.  # m/s      Guestimated cruise speed of Boeing 777
Oswald_factor = 0.76  # Guestimated from ADSEE-I L3

Climb_rate = 2000 * ftmin_to_ms  # m/s      set climb rate by CS25
Climb_gradient = 0.024  # c/V set by CS25 in case 4 engines, should be 0.03
N_engines = 2  # Guestimated

Delta_CD_TO_gear_up = 0.015  # Guestimated from ADSEE-I L3
Delta_CD_TO_gear_down = Delta_CD_TO_gear_up + 0.02  # Guestimated from ADSEE-I L3
Delta_CD_Land = 0.085  # Guestimated from ADSEE-I L3

Delta_oswald_TO = 0.05  # Guestimated from ADSEE-I L3
Delta_oswald_Land = 0.1  # Guestimated from ADSEE-I L3
