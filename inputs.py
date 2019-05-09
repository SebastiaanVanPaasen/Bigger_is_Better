# Put all the inputs that you use in other files in here, the file has been divided in different sections for constants
# input parameters and conversion factors. Always write the units behind a parameter if they exist!
import numpy as np

# Standard constants ---------------------------------------------------------------------------------------------------
Rho_0 = 1.225  # kg/m^3
a = -0.0065  # K/m up to 10 km
g_0 = 9.80565  # m/s^2
T_0 = 288  # K
R_gas = 287  # J/kg/K

# Conversion factors ---------------------------------------------------------------------------------------------------
ftmin_to_ms = 0.00508
kts_to_ms = 0.51444444
lbft2_Nm2 = 47.880172
lbs_to_kg = 0.45359237
ft_to_m = 0.3048
psf_to_nm2 = 1/0.020885
nm_to_km = 1.852
hr_to_sec = 3600
per_hr_to_N = 1 / (hr_to_sec * g_0)

# Define fuel fractions from statistics for the mission profile based on Roskam ----------------------------------------
start = 0.99
taxi = 0.99
t_o = 0.995
climb_1 = 0.98
descent_1 = 0.99
climb_2 = 0.98
descent_2 = 0.99
landing = 0.992

# Passenger characteristics --------------------------------------------------------------------------------------------
N_pas = 450
N_crew = 11
W_person = 205  # lbs

# Class I weight estimation inputs -------------------------------------------------------------------------------------
V_cruise = 499 * kts_to_ms  # m/s  based on the reference aircraft B747
C_fe = 0.003  # estimated from ADSEE-I L3
S_ratio = 6.3  # estimated from ADSEE-I L3
A = 6.96  # based on the reference aircraft B747
Oswald = 0.8  # estimated from ADSEE-I L3
c_j_cruise = 0.75  # 1/hr
c_j_loiter = 0.5  # 1/hr
cruise_range = 1800000  # m based on market analysis
reserve_range = 250 * 1.852 * 1000  # m based on requirement for domestic flights of ADSEE-I L3
W_tfo_frac = 0.003  # estimated from slides ADSEE-I, lecture 3
W_e_frac = 0.525  # Based on average between wide and narrow body, from Ed Obert

# T/W-W/S diagram inputs  ----------------------------------------------------------------------------------------------
Landing_runway = 2500  # m  Guestimated from ADSEE-I L3
Landing_factor = 0.84  # Guestimated from ADSEE-I L3

Rho_TO = Rho_0  # kg/m^3    standard sea-level density
Rho_Cruise = 0.6  # kg/m^3      estimated cruise density
Rho_Landing = Rho_0  # kg/m^3   standard sea-level density

TOP = 220 * lbft2_Nm2  # Guestimated from ADSEE-I L3
Sigma_TO = Rho_TO / Rho_0  # Ratio of densities

V_stall_Cruise = 140  # m/s   Guestimated from ADSEE-I L3
V_stall_Landing = min(np.sqrt(Landing_runway / 0.5847), 65)  # m/s     Guestimated from ADSEE-I L3

CD_0 = 0.018  # guestimated from ADSEE-I L3
CL_Cruise_max = 1.6  # Guestimated from ADSEE-I L3
CL_Landing_max = 3.2  # From Obert
CL_TO_max = 0.8 * CL_Landing_max  # Guestimated from ADSEE-I L3

N_engines = 2  # Guestimated
Climb_rate = 2000 * ftmin_to_ms  # m/s      set climb rate by CS25
Climb_gradient = 0.024  # c/V set by CS25 in case 4 engines, should be 0.03

Delta_CD0_TO_gear_up = 0.015  # Guestimated from ADSEE-I L3
Delta_CD0_TO_gear_down = Delta_CD0_TO_gear_up + 0.02  # Guestimated from ADSEE-I L3
Delta_CD0_Land = 0.085  # Guestimated from ADSEE-I L3

Delta_oswald_TO = 0.05  # Guestimated from ADSEE-I L3
Delta_oswald_Land = 0.1  # Guestimated from ADSEE-I L3
