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
fuel_fractions_input = [start, taxi, t_o, climb_1, descent_1, climb_2, descent_2, landing]

# Define mass fractions from statistics for the aircraft based on Roskam -----------------------------------------------
mass_frac_wing = 0.117
mass_frac_emp = 0.023
mass_frac_fuse = 0.098
mass_frac_nac = 0.018
mass_frac_prop = 0.072
mass_frac_fix = 0.118
mass_fractions_input = [mass_frac_wing, mass_frac_emp, mass_frac_fuse, mass_frac_nac, mass_frac_prop, mass_frac_fix, 0,
                        0, 0]

# Design choices required for the Roskam equations -----------------------------------------------------------------
# first is to include spoilers and speed brakes
# second is with 2 wing mounted engines
# third is with 4 wing mounted engines
# fourth is if landing gear is not wing mounted
# fifth is for strutted wings
# sixth is for fowler flaps
wing_choice = [1, 0, 1, 0, 0, 0]
# choose 1 if you have variable incidence stabilizers
# choose 1 if the horizontal tails are fin mounted
empennage_choice = [0, 0]
# choose 1 if your landing gear is attached to the fuselage
fuselage_choice = 0
# choose 1 if a turbojet or low bypass ratio turbofan is used
nacelle_choice = 1
# choose first if you have buried engine, 0 for no buried engines
# choose secondly if the ducts have a flat cross section, in this case 1
induction_choice = [0, 0]
# choose 1 if piston engines are used
prop_choice = 1
# choose 1 if you have non self-sealing bladder tanks
fuel_sys_choice = 0
# choice is the type of starting system, 1 for one or two jet engines with pneumatic starting system
# 2 for four jet engines with pneumatic starting systems, 3 for jet engines using electric starting systems
# 4 for turboprops with pneumatics, 5 for piston engines using electric systems
start_up_choice = 2
# Choose 1 for fuselage mounted jet, 2 for wing mounted jet, 3 for wing mounted turboprops and 4 for wing mounted
# piston engines the second choice should be 1 if an afterburner is present
engine_choice = [2, 0]
# choice depends on engines, 1 for jet, 2 for turboprop, 3 for radial piston engines, 4 for horizontally opposed
# piston engines
oil_choice = 1
# choice depends on engines, 1 is propellers
hydro_choice = 0

S_ratio = 6.3  # estimated from ADSEE-I L3
C_fe = 0.0045  # estimated from ADSEE-I L3
wing_option = 0  # Depending on the type of wing configuration, 1 is high wing, 0 is low wing
tail_type = 0  # Depending on the type of tail configuration, 1 is T-tail, 0 is conventional

Oswald = 0.9  # estimated from ADSEE-I L3
T_input = 0.27
S_input = 6000.
A = 9
CD_0 = 0.0202
N_engines = 4.  #
w_engine = 3000.  # kg   Obtained from Bram
propeller_choice = 0

W_e_frac_input = 0.525  # Based on average between wide and narrow body, from Ed Obert
CD_cruise_input = (4 / 3) * CD_0
CL_cruise_input = np.sqrt((CD_0 * np.pi * A * Oswald) / 3)

# Inputs that are often changed for the design -------------------------------------------------------------------------
# based on the reference aircraft B747
h_cruise = 8000.  # m based on the sustainability analysis so far
M_cruise = 0.7  # Mach number decided to cruise on

# Passenger characteristics --------------------------------------------------------------------------------------------
N_pas = 450.  # Requirement set by the exercise
N_pas_below = 242.  # Depends on if you want a double-decker 242 if you want a equally sized fuselage for double floor
N_crew = 11.  # Based on the amount of passengers
W_person = 205.  # lbs   Based on statistics of Roskam
W_carg = 20.  # kg Based on statistics

# General aircraft input parameters ------------------------------------------------------------------------------------
# General cruise parameters
h_cruise_input = 8000.  # m based on the sustainability analysis so far
M_cruise_input = 0.7  # Mach number decided to cruise on
Temp_cruise = Temp_0 + a * h_cruise_input  # K  based on the altitude you fly at
a_cruise_input = np.sqrt(gamma * R_gas * Temp_cruise)  # m/s based on the temperature
V_cruise_input = M_cruise_input * a_cruise_input  # m/s  based on the Mach number and speed of sound

# Densities ------------------------------------------------------------------------------------------------------------
Rho_TO = Rho_0  # kg/m^3    standard sea-level density
Rho_Cruise_input = Rho_0 * ((1 + (a * h_cruise_input) / Temp_0) ** (-(g_0 / (R_gas * a))))  # kg/m^3   based on cruise altitude
Rho_Landing = Rho_0  # kg/m^3   standard sea-level density

# aircraft cg-locations ------------------------------------------------------------------------------------------------
x_engines = -1.  # m    x-location engines w.r.t. X_LEMAC
x_fuel = 22.  # m    cg-location fuel w.r.t nose

xcg_oew_mac = 0.25  # m     initial cg location OEW w.r.t. MAC

# engine characteristics -----------------------------------------------------------------------------------------------


# in case propellers are used
prop_characteristics = [2, 50, 1.5, 0.8]  # Number of props, blades per prop, prop diameter and prop efficiency
l_nacelle = 2.  # m   length nacelle

# in case of buried engines
duct_length = 0.
n_inlets = 2.
a_inlets = 3.
n_fuel_tanks = 2.

# main wing ------------------------------------------------------------------------------------------------------------
winglet_height = 0.  # Mostly important for boxed wing, else leave as 0

# horizontal tail ------------------------------------------------------------------------------------------------------
QC_sweep_h = np.radians(32.7)  # degrees    based on wide body statistics
A_h = 4.61  # based on wide body statistics
tap_h = 0.4  # based on statistics from slides ADSEE-I L7
h_tail = np.array([A_h, tap_h, QC_sweep_h])

# vertical tail --------------------------------------------------------------------------------------------------------
QC_sweep_v = np.radians(44.)  # degrees  based on wide body statistics
A_v = 1.73  # based on wide body statistics
tap_v = 0.5  # based on statistics from slides ADSEE-I L7
v_tail = np.array([A_v, tap_v, QC_sweep_v])

# Coefficients ---------------------------------------------------------------------------------------------------------
c_j_cruise = 0.75  # 1/hr
c_j_loiter = 0.5  # 1/hr

# Ranges ---------------------------------------------------------------------------------------------------------------
mission_range = 1400000.  # m     based on market analysis
reserve_range = 250 * 1.852 * 1000  # m based on requirement for domestic flights of ADSEE-I L3
maximum_range = 2000000.  # m    Guestimated

# Fractions ------------------------------------------------------------------------------------------------------------
W_tfo_frac = 0.003  # estimated from slides ADSEE-I, lecture 3

# T/W-W/S diagram inputs  ----------------------------------------------------------------------------------------------
# Landing requirements -------------------------------------------------------------------------------------------------
Landing_runway = 2500.  # m  Guestimated from ADSEE-I L3
Landing_factor = 0.84  # Guestimated from ADSEE-I L3

# Take-off paramters ---------------------------------------------------------------------------------------------------
TOP = 220. * lbft2_Nm2  # Guestimated from ADSEE-I L3
Sigma_TO = Rho_TO / Rho_0  # Ratio of densities

# Stall speeds ---------------------------------------------------------------------------------------------------------
V_stall_Cruise = 250 / 1.2  # m/s   Guestimated from ADSEE-I L3 Take requirements
V_stall_Landing = min(np.sqrt(Landing_runway / 0.5847), 65.)  # m/s     Guestimated from ADSEE-I L3

# Lift and drag coefficients -------------------------------------------------------------------------------------------
CL_Cruise_max = 1.7  # Guestimated from ADSEE-I L3
CL_Landing_max = 3.2  # From Obert
CL_TO_max = 0.9 * CL_Landing_max  # Guestimated from ADSEE-I L3

# Changes in coefficients due to landing or take-off -------------------------------------------------------------------
Delta_CD0_TO_gear_up = 0.015  # Guestimated from ADSEE-I L3
Delta_CD0_TO_gear_down = Delta_CD0_TO_gear_up + 0.02  # Guestimated from ADSEE-I L3
Delta_CD0_Land = 0.085  # Guestimated from ADSEE-I L3

Delta_oswald_TO = 0.05  # Guestimated from ADSEE-I L3
Delta_oswald_Land = 0.1  # Guestimated from ADSEE-I L3

# Estimated parameters -------------------------------------------------------------------------------------------------
Climb_rate = 2000. * ftmin_to_ms  # m/s      set climb rate by CS25
Climb_gradient = 0.024  # c/V set by CS25 in case 4 engines, should be 0.03
