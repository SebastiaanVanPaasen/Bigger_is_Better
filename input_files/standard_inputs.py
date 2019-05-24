import constants_and_conversions as cc

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
mass_fractions = [mass_frac_wing, mass_frac_emp, mass_frac_fuse, 
                        mass_frac_nac, mass_frac_prop, mass_frac_fix, 0,
                        0, 0]

w_tfo = 0.003 # Based on ADSEE
w_e = 0.56 # Based on Roskam

N_pas = 450.  # Requirement
N_pas_below = 450.  # Depends on double-decker config
N_crew = 11.  # Based on the amount of passengers
m_person = 85.  # kg   Based on statistics of Roskam
m_carg = 20.  # kg  Based on statistics

l_nac = 2.  # m     Based on ADSEE
xcg_eng = -1.  # m    Defined w.r.t. LEMAC
xcg_oew_mac = 0.25  #   Defined w.r.t. MAC
x_fuel = 30

N_tanks = 2

CL_cr_max = 1.2  #  Based on ADSEE
