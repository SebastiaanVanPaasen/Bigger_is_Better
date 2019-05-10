# Put all the inputs that you use in other files in here, the file has been divided in different sections for constants
# input parameters and conversion factors. Always write the units behind a parameter if they exist!
from constants_and_conversions import *
import numpy as np

# Define fuel fractions from statistics for the mission profile based on Roskam ----------------------------------------
start = 0.99
taxi = 0.99
t_o = 0.995
climb_1 = 0.98
descent_1 = 0.99
climb_2 = 0.98
descent_2 = 0.99
landing = 0.992
fractions = [start, taxi, t_o, climb_1, descent_1, climb_2, descent_2, landing]

# Passenger characteristics --------------------------------------------------------------------------------------------
N_pas = 450
N_crew = 11
W_person = 205  # lbs

# General aircraft input parameters ------------------------------------------------------------------------------------
# General cruise parameters
h_cruise = 11000  # m based on the sustainability analysis so far
M_cruise = 0.7
T_cruise = T_0 + a * h_cruise
a_cruise = np.sqrt(gamma * R_gas * T_cruise)
V_cruise = M_cruise * a_cruise  # m/s  based on the reference aircraft B747
CL_cruise = 0.7

# Densities ------------------------------------------------------------------------------------------------------------
Rho_TO = Rho_0  # kg/m^3    standard sea-level density
Rho_Cruise = Rho_0 * ((1 + (a * h_cruise) / T_0) ** (-(g_0 / (R_gas * a))))  # kg/m^3   based on cruise altitude
Rho_Landing = Rho_0  # kg/m^3   standard sea-level density

# aircraft cg-locations ------------------------------------------------------------------------------------------------
x_engines = -1.  # m    x-location engines w.r.t. X_LEMAC
x_fuel = 22  # m    cg-location fuel w.r.t nose

xcg_oew_mac = 0.25  # m     initial cg location OEW w.r.t. MAC

# engine characteristics -----------------------------------------------------------------------------------------------
l_nacelle = 2.  # m   length nacelle

# main wing ------------------------------------------------------------------------------------------------------------
airfoil_efficiency = 0.95
alpha = 0.5
alpha_0 = 0.3
CL_des = 1.2
A = 6.96  # based on the reference aircraft B747
S_ratio = 6.3  # estimated from ADSEE-I L3

# horizontal tail ------------------------------------------------------------------------------------------------------
QC_sweep_h = 32.7  # degrees    based on wide body statistics
A_h = 4.61  # based on wide body statistics

# vertical tail --------------------------------------------------------------------------------------------------------
QC_sweep_v = 44  # degrees  based on wide body statistics
A_v = 1.73  # based on wide body statistics

# Coefficients ---------------------------------------------------------------------------------------------------------
C_fe = 0.003  # estimated from ADSEE-I L3
c_j_cruise = 0.75  # 1/hr
c_j_loiter = 0.5  # 1/hr
Oswald = 0.8  # estimated from ADSEE-I L3
CD_0 = C_fe * S_ratio

# Ranges ---------------------------------------------------------------------------------------------------------------
cruise_range = 1800000  # m based on market analysis
reserve_range = 250 * 1.852 * 1000  # m based on requirement for domestic flights of ADSEE-I L3

# Fractions ------------------------------------------------------------------------------------------------------------
W_tfo_frac = 0.003  # estimated from slides ADSEE-I, lecture 3
W_e_frac = 0.525  # Based on average between wide and narrow body, from Ed Obert

# T/W-W/S diagram inputs  ----------------------------------------------------------------------------------------------
# Landing requirements -------------------------------------------------------------------------------------------------
Landing_runway = 2500  # m  Guestimated from ADSEE-I L3
Landing_factor = 0.84  # Guestimated from ADSEE-I L3

# Take-off paramters ---------------------------------------------------------------------------------------------------
TOP = 220 * lbft2_Nm2  # Guestimated from ADSEE-I L3
Sigma_TO = Rho_TO / Rho_0  # Ratio of densities

# Stall speeds ---------------------------------------------------------------------------------------------------------
V_stall_Cruise = 140  # m/s   Guestimated from ADSEE-I L3
V_stall_Landing = min(np.sqrt(Landing_runway / 0.5847), 65)  # m/s     Guestimated from ADSEE-I L3

# Lift and drag coefficients -------------------------------------------------------------------------------------------
CL_Cruise_max = 1.6  # Guestimated from ADSEE-I L3
CL_Landing_max = 3.2  # From Obert
CL_TO_max = 0.8 * CL_Landing_max  # Guestimated from ADSEE-I L3

# Changes in coefficients due to landing or take-off -------------------------------------------------------------------
Delta_CD0_TO_gear_up = 0.015  # Guestimated from ADSEE-I L3
Delta_CD0_TO_gear_down = Delta_CD0_TO_gear_up + 0.02  # Guestimated from ADSEE-I L3
Delta_CD0_Land = 0.085  # Guestimated from ADSEE-I L3

Delta_oswald_TO = 0.05  # Guestimated from ADSEE-I L3
Delta_oswald_Land = 0.1  # Guestimated from ADSEE-I L3

# Estimated parameters -------------------------------------------------------------------------------------------------
N_engines = 2  # Guestimated
Climb_rate = 2000 * ftmin_to_ms  # m/s      set climb rate by CS25
Climb_gradient = 0.024  # c/V set by CS25 in case 4 engines, should be 0.03
