# Add conversion factors or standard constants here if you need them. Make sure to add units behind the constanst!

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

# Constants for control and stability
k_n = -2.5   # -4 = nacelle in front of the LE of the wing or in fuselage nose
             # -2.5 = jet engine pods mounted to de sides of the rear fuselage
V_h_V = 0.85 #  0.85 = fuselage mounted stabilizer 
             #  0.95 = fin mounted stabilizer
             #  1 = for T-tial and canard
eta = 0.95   #  usually used (see slides)
CL_H = -0.8  #lift coefficient for an adjustable tail
             # -1 for full movable tail
             # -0.35*A_H**(1./3.) for a fixed tail