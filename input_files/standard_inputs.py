import constants_and_conversions as cc
import numpy as np

# Ranges of the aircraft
m_range = 1400000.  # m     based on market analysis
r_range = 250 * cc.nm_to_km * 1000  # m     domestic flights of ADSEE-I L3
max_range = 2000000.  # m    Guestimated

# Fuel fraction statistics from Roskam
start = 0.99
taxi = 0.99
t_o = 0.995
climb_1 = 0.98
descent_1 = 0.99
climb_2 = 0.98
descent_2 = 0.99
landing = 0.992
fuel_fractions = [start, taxi, t_o, climb_1, descent_1, climb_2, 
                        descent_2, landing]

# Mass fractions statistics from Roskam 
mass_frac_wing = 0.117
mass_frac_emp = 0.023
mass_frac_fuse = 0.098
mass_frac_nac = 0.018
mass_frac_prop = 0.072
mass_frac_fix = 0.118
mass_frac_lg = 0.02
mass_fractions = [mass_frac_wing, mass_frac_emp, mass_frac_fuse, 
                        mass_frac_nac, mass_frac_prop, mass_frac_fix, mass_frac_lg, 0,
                        0, 0]

w_tfo = 0.003 # Based on ADSEE
w_e = 0.56 # Based on Roskam

N_pas = 450.  # Requirement
N_crew = 11.  # Based on the amount of passengers
m_person = 85.  # kg   Based on statistics of Roskam
m_carg = 20.  # kg  Based on statistics

l_nac = 4.  # m     Based on ADSEE
xcg_eng = -7.  # m    Defined w.r.t. LEMAC
xcg_oew_mac = 0.25  #   Defined w.r.t. MAC
x_fuel = 30  # m    Defined w.r.t. nose

N_tanks = 2
CL_cr_max = 1.2  #  Based on ADSEE

# Stability and Control inputs
V_h_norm = 0.99
s_m = 0.10
SM = 0.10
Cm_0 = -0.1     # airfoil dependent moment coefficient 
CL_0 = 0.5      #  CL of flapped wing at angle of attack of 0 deg 
mu_2 = 1.       # 2D to 3D correction factor from graph  
mu_3 = 0.025    #  correction factor for sweep from graph
dc_c_f = 0.5    # flap geometry ratio, see torenbeek book
x_ac_wing = 0.3

L_run = 2500  # m
V_stall_l = min(np.sqrt(L_run / 0.5847), 65.)
CL_l_max = 3.2 

N_cargo = 2     # number of cargo compartments
cargo_fwdfrac = 0.6     # Fraction of the amount the front compartment holds, if N_cargo = 1 then the value is 1
delta_flap =  0.5235    # flap deflection [rad]
alpha_land =   0.1745   # angle of attach during approach [rad]
CL_H        = -0.8  # liftcoefficient adjustable tail
eta        =  0.95  # airfoil efficiency = 0.95
x_ac       =  0.3   # aerodynamic centre location over the mac [n/m]

