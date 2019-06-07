import numpy as np
import constants_and_conversions as cc

from loading_definitions import Section, Helpers
from scipy.integrate import quad
from read_csv_input import read_output

filename = 'HIGH SEMIDD 2E STRUT'
weights, wing, cruise_conditions = read_output(filename)


S = wing["S"]
b = wing["b"]
AR = wing["A"]
taper = wing["Taper"]
qcsweep= wing["Sweep"]
Cr = wing["C_root"]
Ct = wing["C_tip"]
 
W_TO = weights["W_TO"]
W_W = weights["W_W"]
M_F = (weights["W_F"] / cc.g_0) / 2
N_eng = 2
W_eng = weights["W_E"] + weights["W_N"] / N_eng


T_TO = cruise_conditions["T_TO"]
alt = cruise_conditions["H_cr" ]
V_cr = cruise_conditions["V_cr"] * 1.3
CD0 = cruise_conditions["CD_0"]


sweep_LE = np.arctan(np.tan(qcsweep) - 4 / AR * (-0.25 * (1 - taper) / (1 + taper))) # rad
sweep_SC = np.arctan(np.tan(sweep_LE) - 4 / AR * (0.4 * (1 - taper) / (1 + taper)))  # sweep shear center
tc = 0.14
tot_volume = (0.5 * ((Cr ** 2) * tc + (Ct ** 2) * tc) * (b / 2)) * 2  # m^3
w_spec_weight = W_W / tot_volume  # N/m^3
w_spec_fuel = 0.804 * 1000 # kg/m^3



pos_thrust = -1.5  # position of y of the thrust vector
pos_eng_1 = 7.
#pos_eng_2 = 0.6 * (b / 2)
pos_eng = [pos_eng_1]
pos_strut = -1
strutforce = 1000000.
lift_pos, weight_pos = 1/4, 1/2
fuel_pos, strut_pos, eng_pos = 1/2, 1/4, 0

force_applications = [lift_pos, weight_pos, fuel_pos, eng_pos, 0.5]

Rho_cr = cc.Rho_0 * ((1 + (cc.a * alt) / cc.Temp_0) ** (-(cc.g_0 / (cc.R_gas * cc.a))))
n_ult= 3.75
n = n_ult / 1.5


x = 0
start_fuel = 15.
x_strut = 14.
d_fus = 6.8
strut_applic = 0.25
pos_strut, l_strut, T_y, T_z = Helpers.calc_strut_char(Cr, Ct, b, d_fus, x_strut, 
                                             sweep_LE, strut_applic, tc)
E = 69 * (10 ** 9)
E_strut = 181 * (10 ** 9)
A_strut = (1 / 4) * np.pi * ((10 / 100) ** 2)
print(l_strut)

dx = 0.2
discr = np.arange(0, b / 2, dx) + dx / 2
cl = Helpers.input_CL(S, V_cr, Rho_cr, W_TO)
polar_cl, polar_cd = Helpers.calc_polars(cl)

moment_forces = []
num = []
lift_force = []
weight_force = []

T = np.zeros((1, len(discr)))
F_y = np.zeros((1, len(discr)))
F_z = np.zeros((1, len(discr)))
M_z = np.zeros((1, len(discr)))
M_y = np.zeros((1, len(discr)))
    
for i in range(len(discr)):
#    print("Starting on position ", discr[i])
   
    section = Section(Cr, Ct, b, discr[i], dx, tc, n, CD0, pos_eng_1)
    section.run_geometrics(sweep_LE)
    
    lift, weight, drag = section.calc_forces(w_spec_weight, polar_cl, polar_cd, Rho_cr, V_cr)
    M_F, fuel_weight = section.calc_fuel_weight(w_spec_fuel, M_F, start_fuel)
    
    section.calc_engine_char(pos_eng, T_TO, W_eng, N_eng)
    W_engine = section.calc_torques(force_applications, sweep_LE, sweep_SC, pos_thrust)

    T[0][i] = section.calc_torques(force_applications, sweep_LE, sweep_SC, pos_thrust)
    F_y[0][i], F_z[0][i] = section.calc_tot_force()  
    M_z[0][i], M_y[0][i] = section.calc_moments()
    
    
    moment_forces.append(lift)
    lift_force.append(lift)
    weight_force.append(weight)
    num.append(-lift * discr[i])
    
    moment_forces.append(-weight)
    num.append(weight * discr[i])
    
    moment_forces.append(-W_engine)
    num.append(W_engine * pos_eng_1)
#
#
print((lift_force[0] + lift_force[-1]) / 2)
print((weight_force[0] + weight_force[-1]) / 2)
print(W_eng)
print(lift_force[0], lift_force[-1])
print(weight_force[0], weight_force[-1])
#T = np.sum(T)   
#F_y = np.sum(F_y) 
#F_z = np.sum(F_z) 
#M_y = np.sum(M_y) 
#M_z = np.sum(M_z) 
#
#M_forces = np.sum(np.array(moment_forces))
#M_num = np.sum(np.array(num))
#
#M_forces = M_forces
#M_num = M_num
#
#print(T, F_y, F_z, M_y, M_z, M_forces, M_num)
#
#forces = [T, F_y, F_z, M_y, M_z, M_forces, M_num]
#
#I = 0
#i = 0
#while discr[i] < x_strut:
#    chord = Helpers.calc_chord(Cr, Ct, b, discr[i])
#    part_inertia = Helpers.calc_inertia(chord * 0.1, (1 / 100), 0.5 * tc * chord)
#    I += part_inertia
#    i += 1
#
#print(I / len(discr))
#
##sigma_strut = 552 * (10 ** 6)
#force_results = Helpers.indet_sys(pos_strut, E_strut, l_strut, A_strut, E, I, x_strut, forces, T_y, T_z)
#
#print(force_results)
#
##print((sigma_strut * l_strut) / (E_strut) * 1000)
#print((force_results[0] * l_strut) / (E_strut * A_strut) * 1000)
#
#print((force_results[0] / A_strut) / (10 ** 6))
#
#m_z = force_results[6]
#v_y = force_results[2]




#sum_x = np.array([pos_strut[0], 1, 0, 0])
#sum_y = np.array([pos_strut[1], 0, 1, 0])
#sum_m = np.array([-2*x_strut, 0, 0, 1])




    




#T = np.zeros((1, len(discr)))
#F_y = np.zeros((1, len(discr)))
#F_z = np.zeros((1, len(discr)))
#M_z = np.zeros((1, len(discr)))
#M_y = np.zeros((1, len(discr)))
#D = []
#
#
#for i in range(len(discr)):
#    print("Starting on position ", discr[i])
#    
#    section = Section(Cr, Ct, b, discr[i], dx, tc, n, CD0, pos_eng_1)
#    
#    section.run_geometrics(sweep_LE)
#    drag = section.calc_forces(w_spec_weight, polar_cl, polar_cd, Rho_cr, V_cr)
#    
#    D.append(drag)
#    
#    M_F = section.calc_fuel_weight(w_spec_fuel, M_F)
#    
#    section.calc_strut_force(pos_strut, strutforce)
#    section.calc_engine_char(pos_eng, T_TO, W_eng, N_eng)
#          
#    T[0][i] = section.calc_torques(force_applications, sweep_LE, sweep_SC, pos_thrust)
#    F_y[0][i], F_z[0][i] = section.calc_tot_force()  
#    M_z[0][i], M_y[0][i] = section.calc_moments()
#    
#    x += dx
#
#   
#discr = np.append(np.array([0]), discr)
#
#T_root = np.sum(T)   
#F_y_root = np.sum(F_y) 
#F_z_root = np.sum(F_z) 
#M_y_root = np.sum(M_y) 
#M_z_root = np.sum(M_z) 
#
#
#T_plot = [T_root]
#F_y_plot = [F_y_root]
#F_z_plot = [F_z_root]
#M_z_plot = [M_z_root]
#M_y_plot = [M_y_root]
#
#    
#for i in range(len(discr) - 1):
#    T_plot.append(T_plot[i]-T[0][i])
#    F_y_plot.append(F_y_plot[i] - F_y[0][i])
#    F_z_plot.append(F_z_plot[i] - F_z[0][i])
#    
#    
#for i in range(len(discr) - 1):
#    M_y_plot.append(M_y_plot[i] + F_z_plot[i + 1] * (discr[i + 1] - discr[i]))
#    M_z_plot.append(M_z_plot[i] - F_y_plot[i + 1] * (discr[i + 1] - discr[i]))
#
#Helpers.plotter(discr, T_plot, M_z_plot, M_y_plot, F_y_plot, F_z_plot)
#
##print(T_root)
##print(F_y_root)
##print(F_z_root)
##print(M_y_root)
##print(M_z_root)
#
##print(F_y_plot)
##print(M_y)    
##print(M_y_plot)  
