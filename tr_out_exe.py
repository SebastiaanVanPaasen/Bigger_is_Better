import numpy as np
import constants_and_conversions as cc
import matplotlib.pyplot as plt

from try_out import Section, Helpers
from read_csv_input import read_output

filename = 'HIGH SEMIDD 2E STRUT'
weights, wing, cruise_conditions = read_output(filename)


CD0 = cruise_conditions["CD_0"]
S = wing["S"]
b = wing["b"]
AR = wing["A"]
taper = wing["Taper"]
qcsweep= wing["Sweep"]
Cr = wing["C_root"]
Ct = wing["C_tip"]
 
W_TO = weights["W_TO"]
W_W = weights["W_W"]
M_F = weights["W_F"] / cc.g_0
N_eng = 2
W_eng = weights["W_E"]+ weights["W_N"] / N_eng


T_TO = cruise_conditions["T_TO"]
alt = cruise_conditions["H_cr" ]
V_cr = cruise_conditions["V_cr"] * 1.3


sweep_LE = np.arctan(np.tan(qcsweep) - 4 / AR * (-0.25 * (1 - taper) / (1 + taper))) # rad
sweep_SC = np.arctan(np.tan(sweep_LE) - 4 / AR * (0.4 * (1 - taper) / (1 + taper)))  # sweep shear center
tc = 0.14
tot_volume = (0.5 * ((Cr ** 2) * tc + (Ct ** 2) * tc) * (b / 2)) * 2  # m^3
w_spec_weight = W_W / tot_volume  # N/m^3
w_spec_fuel = 0.804 * 1000 # kg/m^3



pos_thrust = -1.5  # position of y of the thrust vector
pos_eng_1 = 7
#pos_eng_2 = 0.6 * (b / 2)
pos_eng = [pos_eng_1]
pos_strut = 15
strutforce = 1000000.
lift_pos, weight_pos, fuel_pos, eng_pos = 1/4, 1/2, 1/2, 0

force_applications = [lift_pos, weight_pos, fuel_pos, lift_pos, eng_pos]

Rho_cr = cc.Rho_0 * ((1 + (cc.a * alt) / cc.Temp_0) ** (-(cc.g_0 / (cc.R_gas * cc.a))))
n_ult= 3.75
n = n_ult / 1.5


x = 0
start = 0
dx = 1
discr = np.arange(0, b / 2, dx)
cl = Helpers.input_CL(S, V_cr, Rho_cr, W_TO)
polar_cl, polar_cd = Helpers.calc_polars(cl)

T = np.zeros((1, len(discr)))
F_y = np.zeros((1, len(discr)))
F_z = np.zeros((1, len(discr)))
M_z = np.zeros((1, len(discr)))
M_y = np.zeros((1, len(discr)))

for i in range(len(discr)):
    print("Starting on position ", discr[i])
    
    section = Section(Cr, Ct, b, discr[i], dx, tc, n, CD0)
    
    section.run_geometrics(sweep_LE)
    section.calc_forces(w_spec_weight, polar_cl, polar_cd, Rho_cr, V_cr)
    
    M_F = section.calc_fuel_weight(w_spec_fuel, M_F)
    
    section.calc_strut_force(pos_strut, strutforce)
    section.calc_engine_char(pos_eng, T_TO, W_eng, N_eng)
    
    T[0][i] = section.calc_torques(force_applications, sweep_LE, sweep_SC, pos_thrust)
    F_y[0][i], F_z[0][i] = section.calc_tot_force()
    M_z[0][i], M_y[0][i] = section.calc_moments(start, force_applications)
    
    x += dx
    
print(T)
print(F_y)
print(F_z)
print(M_z)
print(M_y)

plt.subplot(2,3,5)
plt.subplot(2,3,1)
plt.gca().set_title('Mz distribution')
plt.plot(discr, M_z[0])
plt.xlabel('Position Along Wing Span [$m$]')
plt.ylabel('Mz $[Nm]$')
plt.subplot(2,3,2)
plt.gca().set_title('Fz distribution')
plt.plot(discr, F_z[0])
plt.xlabel('Position Along Wing Span [$m$]')
plt.ylabel('Fz [$N$]')
plt.subplot(2,3,3)
plt.gca().set_title('My distribution')
plt.plot(discr, M_y[0])
plt.xlabel('Position Along Wing Span [$m$]')
plt.ylabel('My $[Nm]$')
plt.subplot(2,3,4)
plt.gca().set_title('Fy distribution')
plt.plot(discr, F_y[0])
plt.xlabel('Position Along Wing Span [$m$]')
plt.ylabel('Fy [$N$]')
plt.subplot(2,3,5)
plt.gca().set_title('T distribution')
plt.plot(discr, T[0])
plt.xlabel('Position Along Wing Span [$m$]')
plt.ylabel('T $[Nm]$')
plt.show()