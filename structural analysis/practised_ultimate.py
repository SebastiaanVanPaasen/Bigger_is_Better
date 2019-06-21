# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 11:13:33 2019

@author: mathi
"""

import numpy as np
import matplotlib.pyplot as plt
import constants_and_conversions as cc
import centroid_wing as cw
import Airfoil_inertia as ai
import parameter_requirements as pr

from scipy.interpolate import interp1d
from class_I.lift_distr import get_correct_data, lift_distribution

AR = 15
D_fus = 7.3

#R_strut = 5 / 1000
#A_strut = 0.25 * np.pi * ((2 * R_strut) ** 2) ##
E_strut = 60.1 * (10 ** 9)  
AR = 15
taper = 0.357 
MAC = 4.03

cr = 5.53
ct = 1.97
b = 56.3
S = 211.19
dihedral = (1.5/180) * np.pi

LE_root = 19.11
pos_of_chord = 0.25
sweep = 0
pos_from_centerline = 2.5
strut_loc_fus = LE_root + pos_of_chord * cr + np.tan(sweep)* (b/2)
strut_heigth = D_fus - np.tan(dihedral)*(b/2) - 1.

W_wing = 156375
E_wing = 75 * (10 ** 9)
I_zz_wing = 0.25
L_wing = b / 2

H_cr = 9000
V_cr = 218
rho_cr = cc.Rho_0 * ((1 + (cc.a * H_cr) / cc.Temp_0) ** (-(cc.g_0 / (cc.R_gas * cc.a) + 1)))
CD_0_cr = 0.02114
print(rho_cr)
W_TO = 1520276 
W_fuel = 170698
W_N = 25875

N_eng = 2
W_eng_v = 102023/2 #4760 * cc.g_0 #+ W_N / N_eng
Rho_fuel = 0.804 * 1000


def deter_d_force(applic, x, force, a_s, I_wing):
    
#    I_wing = I_wing[::-1]
    Y_force = (force/ (6 * E_wing * I_wing)) * (2 * ((L_wing - a_s) ** 3) - 3 * ((L_wing - a_s) ** 2) * applic + applic ** 3)
    Theta_force= (-force / (2 * E_wing * I_wing)) * (((L_wing - a_s) - applic) ** 2) * x
    
    Mac = np.array([])
#    print("doei",I_wing)
    for i in range(len(Theta_force)):
        
        
        if x[i] <= applic:
            Mac = np.append(Mac, 0)
                
        else:
            mac = (force / (6 * E_wing * I_wing[i])) * ((x[i] - applic) ** 3)
            Mac = np.append(Mac, mac)
    

    d_force = Y_force + Theta_force + Mac
    
#    plt.plot(x, Y_force)
#    plt.plot(x, Theta_force)
#    plt.plot(x, Mac)
#    plt.show()
    
    return d_force


def get_data():
#    make_avl_file()
    
    output_avl = lift_distribution()
    x_pos, cl, cdi = get_correct_data(output_avl, MAC)
#    print(x_pos, cl, cdi)

    PolyFitCurveCl = interp1d(x_pos, cl, kind="cubic", fill_value="extrapolate")
    PolyFitCurveidrag = interp1d(x_pos, cdi, kind='cubic', fill_value='extrapolate')
    
    return PolyFitCurveCl,  PolyFitCurveidrag


def calc_chord(x):
    return cr - ((cr - ct) / L_wing) * x
    
#print(ai.s_airfoil(100,60, calc_chord))
    
def calc_area(x, width):
    chord = calc_chord(x)
    
    return width * chord
        

def deter_lift(cl_curve, x, width):
    CL = np.array([])
    S = np.array([])
    
    for i in range(len(x)):
        CL = np.append(CL, cl_curve(x[i]))
        S = np.append(S, calc_area(x[i], width))
        
      
    qs = 0.5 * rho_cr * (V_cr ** 2) * S
    L = qs * CL
    
#    plt.plot(x, L)
#    plt.show()
    
    return L


def deter_drag(cd_curve, x, width, n):
    CD_nom = np.array([])
    CD = np.array([])
    S = np.array([])
    
    for i in range(len(x)):
        CD_nom = np.append(CD_nom, cd_curve(x[i]) + CD_0_cr)
        CD = np.append(CD, n * cd_curve(x[i]) + CD_0_cr)
        S = np.append(S, calc_area(x[i], width))
        
        
    qs = 0.5 * rho_cr * (V_cr ** 2) * S
    qs_nom = 0.5 * rho_cr * (V_cr ** 2) * S
    
    D = qs * CD
    D_nom = qs_nom * CD_nom
    
#    plt.plot(x, D)
#    plt.plot(x, D_nom)
#    plt.show()
    
    return D, np.sum(D_nom)


def deter_weight(w_wing, x, width):
    W = np.array([])
    Volumes = np.array([])
    
    for i in range(len(x)):
        Volumes = np.append(Volumes, (0.6 * 0.0809 * (calc_chord(x[i]) ** 2) * width))
        
        
    V_tot = np.sum(Volumes)
    w_spec = w_wing / V_tot
    

    for i in range(len(x)):
        W = np.append(W, Volumes[i] * w_spec * 1.02)
        
#    plt.plot(x, W)
#    plt.show()
    
    return W, Volumes


def deter_fuel(w_fuel, volumes, density, x, x_start):
    W = np.array([])
#    print("total volume", sum(volumes))
#    print(w_fuel)
    for i in range(len(x)):
        
        
        if x[i] > x_start and w_fuel > 0:
            w_fuel_section = density * volumes[i] * cc.g_0 
        else:
            w_fuel_section = 0


        W = np.append(W, w_fuel_section)
        w_fuel -= w_fuel_section
    
#    print(W)
#    print(sum(W))
#    plt.figure()
#    plt.plot(x, W, label = "start fuel tank " + str(round(x_start,2)))
#    plt.xlabel("X-position [m]")
#    plt.ylabel("Fuel weigth [N]")
#    plt.show()
#    plt.legend()
    
    return W
            
    
def indet_sys(dx, angle, L_s, a_s, a_e, cl_polar, I_zz_sections):   
    X_root = np.arange(0, L_wing+dx, dx)
 
    X_tip = np.arange(L_wing - dx / 2, 0, -dx)


    x_start = L_wing * 0.28
    
    alpha = np.radians(0.9153)
    n = -1*1.5
    
    
    lifts_v = n * deter_lift(cl_polar, X_tip[::-1], dx)   
#    print("weight", W_TO)
#    print("total lift", 2*sum(lifts_v))
    drags_v, thrust = deter_drag(cd_polar, X_tip[::-1], dx, n)   
    weights_v, Volumes = deter_weight(W_wing, X_tip[::-1], dx)  #
    fuel_weights_v = deter_fuel(W_fuel / 2, Volumes, Rho_fuel, X_tip[::-1], x_start)
    
    lifts = np.cos(alpha) * lifts_v - np.sin(alpha) * drags_v
    drags = np.sin(alpha) * (-lifts_v + weights_v + fuel_weights_v + W_eng_v) + np.cos(alpha) * drags_v 
    
    weights = np.cos(alpha) *  weights_v
    fuel_weights = np.cos(alpha) * fuel_weights_v
    W_eng = np.cos(alpha) * W_eng_v
    

#    print("Forces of lift, drag, thrust, weight and fuel weight")
#    print(np.sum(lifts))
#    print(np.sum(drags))
#    print(np.sum(thrust))
#    print(np.sum(weights))
#    print(np.sum(fuel_weights))
#    print()
    
    defl_l = np.zeros(len(X_root))
    defl_w = np.zeros(len(X_root))
    defl_wf = np.zeros(len(X_root))
    defl_eng = np.zeros(len(X_root))

    
    for i in range(len(X_tip)):
        defl_l += deter_d_force(X_tip[i], X_root, lifts[i], 0, I_zz_sections[::-1])
        defl_w += deter_d_force(X_tip[i], X_root, -weights[i], 0, I_zz_sections[::-1])
        defl_wf += deter_d_force(X_tip[i], X_root, -fuel_weights[i], 0, I_zz_sections[::-1])
        
      
    defl_eng += deter_d_force(A_E, X_root, -W_eng, 0, I_zz_sections[::-1])
    

    sigma = 140 * (10 ** 6)
    d_wing = defl_l + defl_w + defl_wf + defl_eng
    
#    print("d wing")
#    print(d_wing[int(a_s / dx)])
#    print("lift", defl_l[int(a_s / dx)])
#    print("weight", defl_w[int(a_s / dx)])
#    print("fuel", defl_wf[int(a_s / dx)])
#    print("eng", defl_eng[int(a_s / dx)])
#    print()
    
#    print("strut length", L_s)
#    print("angle of strut", np.sin(angle), angle)
#    print("sigma", sigma)
    d_strut = np.sin(angle) * sigma * L_s / E_strut 
    print("deflection of strut", d_strut)
    d_strut_v = d_wing[int((a_s) / dx)] - d_strut
    print(d_strut_v)
    y_a = -(2 * L_wing ** 3 - 3 * (L_wing ** 2) *  (a_s + (dx/2)) +  (a_s + (dx/2)) ** 3)/(6 * E_wing * I_zz_sections[::-1][int((a_s)/dx)])
    theta_a = (((L_wing - (a_s + dx/2)) ** 2) / (2 * E_wing * I_zz_sections[::-1][int((a_s)/dx)])) * (a_s+(dx/2))
    F_strut = -(d_strut_v / (y_a + theta_a))  # (2 * ((L_wing) ** 3))) * (6 * E_wing * I_wing)) / np.sin(angle)
    distributions = [lifts, weights, fuel_weights, W_eng, drags, thrust]
    
    return F_strut, d_strut, d_strut_v, d_wing, distributions

    
def strut_opt(A_S, A_E, cl_curve, width, I_wing, gamma, L_strut):
    
#    print("hoi",I_wing[int((A_S) / width)])
    F_strut = np.arange(0, 6000000, 1000)
    force, deflection, all_forces = indet_sys(F_strut, width, gamma, L_strut, A_S, A_E, cl_curve, I_wing[int((A_S) / width)])
    
#        print("First found optimum")
#        print(force, deflection)
#        print()
    
    
    F_strut = np.arange(force - 2000, force + 2000, 0.1)
    strut_force, deflection, all_forces = indet_sys(F_strut, width, gamma, L_strut, A_S, A_E, cl_curve, I_wing[int((A_S) / width)])
    
#    print(deflection)
#    print("Final optimum")
#    print(force, deflection)
#    print()
   
    return strut_force, deflection, all_forces


A_E = 21.16 #7m from center of fuselage
A_S_L = np.arange(8.05, 8.50, 0.5)

cl = W_TO / (0.5 * rho_cr * (V_cr ** 2) * S)
#print(cl)
cl_polar, cd_polar = get_data()

dx = 0.1

X_root = np.arange(0, L_wing+dx, dx)

X_tip = np.arange(L_wing - dx / 2, 0, -dx)

#I_zz_sections = np.zeros(len(X_root))
Mz_dist = np.zeros((len(A_S_L), len(X_root)))
My_dist = np.zeros((len(A_S_L), len(X_root)))
Vy_dist = np.zeros((len(A_S_L), len(X_root)))
Vz_dist = np.zeros((len(A_S_L), len(X_root)))

print("-------Starting on strut optimisation-------")

N = (L_wing / dx) 
l_spar_h, t_spar_v, t_spar_h = cw.l_spar_h, cw.t_spar_v, cw.t_spar_h
boom_area_old = 0.005
boom_area_all = np.zeros(len(A_S_L)) 
 
F_strut = np.zeros(len(A_S_L))
L_str = np.zeros(len(A_S_L))

for idx in range(len(A_S_L)):

#    I_zz_spar, I_yy_spar, I_yz_spar = ai.I_zz_spars(l_spar_h, t_spar_v, t_spar_h, N, b ,calc_chord, boom_area_old,X_root, dx)
#    I_zz_req = pr.required_Izz(N, b, calc_chord, Mz_dist[idx], boom_area_old, X_root, dx)
#
#    boom_area_new = ai.wing_geometry(I_zz_req, I_zz_spar, N, b, calc_chord, boom_area_old, X_root, dx)

#    print(abs(boom_area_new - boom_area_old) )
#    while abs(boom_area_new - boom_area_old) > 1 / 100000:
#        print(abs(boom_area_new - boom_area_old) )
#        print("New iteration")
#    boom_area_old = boom_area_new
    boom_area_all = 0.0045
    
    I_zz_spar, I_yy_spar, I_yz_spar = ai.I_zz_spars(l_spar_h, t_spar_v, t_spar_h, N, b ,calc_chord, boom_area_all,X_root, dx)
    I_zz_req = pr.required_Izz(N, b, calc_chord, Mz_dist[idx], boom_area_all, X_root, dx)
    
#        boom_area_new = ai.wing_geometry(I_zz_req, I_zz_spar, N, b, calc_chord, boom_area_old,X_root, dx)
#        airfoil_area, z_c_airfoil, y_c_airfoil = cw.get_skin_centroid(N, b, calc_chord,dx)
#
#        z_centroid_all_sec, y_centroid_all_sec, y_loc_spar_up, y_loc_spar_low, y_loc_stiff_up, y_loc_stiff_low, y_vertical_spar, z_loc_stiff_up, spar_loc_sec, z_loc_stiff_low, spar_areas_verti = cw.wing_centroid(boom_area_new, cw.spar_areas_hori, cw.t_spar_v, z_c_airfoil, y_c_airfoil, cw.n_stiff_up, cw.n_stiff_low, N, b, calc_chord, X_root, dx)
#
#    
#    #    print("I_zz_req",I_zz_req)
#        boom_area_all[idx] = boom_area_new[0]
#        print("Updated boom area for strut pos" + str(A_S_L[idx]))
#        
#        print(boom_area_all[idx] * 10000)
#        print()
#    
    
    
    I_zz_sections, I_yy_wing, I_yz_wing = ai.inertia_wing(I_zz_spar, I_yy_spar, I_yz_spar, boom_area_all, N, b, calc_chord, X_root, dx)
#       
    print("I_zz_sections",I_zz_sections)
#        I_zz_sections = np.array(293 * [0.00499])
    strut_loc_wing = L_wing - A_S_L[idx] - pos_from_centerline
   
    gamma = np.arctan(strut_heigth / strut_loc_wing)
    L_strut = np.sqrt(strut_heigth**2 + strut_loc_wing**2)
    L_str[idx] = L_strut
    print("length strut", L_strut)


    F_str, d_str, d_str_v, d_w, all_forces = indet_sys(dx, gamma, L_strut, A_S_L[idx], A_E, cl_polar, I_zz_sections)
#    print("Strut force and deflections")
#    print("Strut force", F_str)
#    print("Deflection of the strut", d_str)
##        print(d_str_v)
#    print("Required strut area", F_str / (100 * (10 ** 6)))
##        print()
#    F_strut[idx] = F_str
##    F_str = results[0]
#    print("Strut force in first optimisation")
#    print(F_str)
#    print()
    
    Lift, Weight, Fuel_weight, W_eng, Drag, Thrust = all_forces

    Lift_mom = Lift * X_tip[::-1]
    Weight_mom = Weight * X_tip[::-1]
    Fuel_mom = Fuel_weight * X_tip[::-1]
    Eng_mom = W_eng * (L_wing - A_E)
#    Strut_mom = F_str * (L_wing - A_S_L[idx])
#        Strut_mom_hori = F_str/np.tan(angle) * y_centroid_all_sec[int((L_wing - A_S_L[idx])/dx)] - 
    Drag_mom = Drag * X_tip[::-1]
    Thrust_mom = Thrust * (L_wing - A_E)
    
    Mz_root = sum(Lift_mom) - sum(Weight_mom) - sum(Fuel_mom) - Eng_mom #- Strut_mom
    My_root = sum(Drag_mom) - Thrust_mom
    Vy_root = sum(Lift) - sum(Weight) - sum(Fuel_weight) - W_eng #- F_str
    Vz_root = Thrust - sum(Drag)
    
    print("Root moment around z: ", Mz_root)
    print("Root moment around y: ", My_root)
    print("Root force around y: ", Vy_root)
    print("Root force around z: ", Vz_root)
    print()
#        
    Mz_dist[idx][0] = Mz_root
    My_dist[idx][0] = My_root
    Vy_dist[idx][0] = Vy_root
    Vz_dist[idx][0] = Vz_root
    
    for i in range(len(X_tip)):
        Vy_section = Vy_dist[idx][i] - Lift[i] + Weight[i] + Fuel_weight[i]
        Vz_section = Vz_dist[idx][i] + Drag[i]
        
        if X_root[i] > (L_wing - A_E) and X_root[i - 1] < (L_wing - A_E):
            Vy_section += W_eng
            
        if X_root[i] > (L_wing - A_S_L[idx]) and X_root[i - 1] < (L_wing - A_S_L[idx]):
            Vy_section += F_str
            
        if X_root[i] > (L_wing - A_E) and X_root[i - 1] < (L_wing - A_E):
            Vz_section -= Thrust
        
        Vy_dist[idx][i + 1] = Vy_section
        Vz_dist[idx][i + 1] = Vz_section
    
    for i in range(len(X_tip)):
        Mz_dist[idx][i + 1] = Mz_root - Vy_root * X_root[i+1] + sum(((Lift - Weight - Fuel_weight) * (X_root[i+1] - X_tip[::-1]))[:i])
        My_dist[idx][i + 1] = My_root + Vz_root * X_root[i+1] + sum((Drag*(X_root[i+1] - X_tip[::-1]))[:i])
        
        if X_root[i+1] > (L_wing - A_E): #and X_root[i - 1] < (L_wing - A_E):
            Mz_dist[idx][i + 1] -= W_eng * (X_root[i+1]- (L_wing - A_E))
            
        if X_root[i+1] > (L_wing - A_S_L[idx]): #and X_root[i - 1] < (L_wing - A_S_L[idx]):
            Mz_dist[idx][i + 1] -= F_str *(X_root[i+1]- (L_wing - A_S_L[idx]))
            
        if X_root[i+1] > (L_wing - A_E): #and X_root[i - 1] < (L_wing - A_E):

            My_dist[idx][i + 1] -= Thrust * (X_root[i+1]- (L_wing - A_E))
    d_lift = 0
    d_weight = 0
    d_fuel_weight = 0
    diff = []
    diff_2 = []
    
    for i in range(len(X_tip)):
        d_lift += deter_d_force(X_tip[i], X_root, Lift[i], 0, I_zz_sections[::-1])
        d_weight += deter_d_force(X_tip[i], X_root, -Weight[i], 0, I_zz_sections[::-1])
        d_fuel_weight += deter_d_force(X_tip[i], X_root, -Fuel_weight[i], 0, I_zz_sections[::-1])
        

    d_strut = deter_d_force(A_S_L[idx] + dx/2, X_root, -F_str, 0, I_zz_sections[::-1])
    d_engine = deter_d_force(A_E, X_root, -W_eng, 0, I_zz_sections[::-1])
    

    d = d_lift + d_weight + d_fuel_weight + d_strut + d_engine
    
#    x_start = 0*(b/2)
#    plt.figure()
#    plt.plot(X_tip[::-1], all_forces[2], label = "start fuel tank " + str(round(x_start,2)))
#    plt.xlabel("X-position [m]")
#    plt.ylabel("Fuel weigth [N]")
#    plt.show()
#    plt.legend()
#        
#        for i in range(1, len(d)):
#            diff.append(d[i - 1] - d[i])
#            
# 
#        for i in range(1, len(diff)):
#            diff_2.append(diff[i - 1] - diff[i])
            

#    print("deflection of components after strut calculation")
#    print(d_lift[int((A_S_L[idx])/dx)])
#    print(d_weight[int((A_S_L[idx])/dx)])
#    print(d_fuel_weight[int((A_S_L[idx])/dx)])
#    print(d_engine[int((A_S_L[idx])/dx)])
#    print(d_strut[int((A_S_L[idx])/dx)])
#    print(d[int((A_S_L[idx])/dx)])
#    print("izz", I_zz_sections[::-1])
    print(Vy_dist[idx][0])
    print(Mz_dist[idx][0])
    print(My_dist[idx][0])
    print(Vz_dist[idx][0])
#    
#    plt.figure()
    plt.subplot(2, 3, 1)
    plt.plot(X_root, Vy_dist[idx], label = "Vy for pos " + str(A_S_L[idx]))
    plt.xlabel("X-position [m]")
    plt.ylabel("Vy [N]")
    plt.title("Vy distribution")
    ax = plt.gca()
    ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
#    ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
#    plt.legend()
    
    plt.subplot(2, 3, 4)
    plt.plot(X_root, Mz_dist[idx], label = "Mz for pos " + str(A_S_L[idx]))
    plt.xlabel("X-position [m]")
    plt.ylabel("Mz [Nm]")
    plt.title("Mz distribution")
    ax = plt.gca()
    ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))

#    plt.legend()
    
    plt.subplot(2, 3, 2)
    plt.plot(X_root, Vz_dist[idx], label = "Vz for pos " + str(A_S_L[idx]))
    plt.xlabel("X-position [m]")
    plt.ylabel("Vz [N]")
    plt.title("Vz distribution")
    ax = plt.gca()
    ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))

#    plt.legend()
    
    plt.subplot(2, 3, 3)
#    plt.plot(X_tip, d_lift, label = "lift for pos " + str(A_S_L[idx]))
#    plt.plot(X_tip, d_weight, label = "weight for pos " + str(A_S_L[idx]))
#    plt.plot(X_tip, d_fuel_weight, label = "fuel weight for pos " + str(A_S_L[idx]))
#    plt.plot(X_tip, d_strut, label = "strut for pos " + str(A_S_L[idx]))
#    plt.plot(X_tip, d_engine, label = "engine for pos " + str(A_S_L[idx]))
    plt.plot(X_root[::-1], d, label = "Deflection for pos " + str(A_S_L[idx]))
    plt.xlabel("X-position [m]")
    plt.ylabel("Deflection [m]")
#    plt.title("Deflection of loading components along the span")
    plt.title("Deflection along the span")
#    ax = plt.gca()
#    ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))

#    plt.legend()
#    plt.plot(X_root[::-1], d_lift, label = "lift ")#for pos " + str(A_S_L[idx]))
#    plt.plot(X_root[::-1], d_weight, label = "weight")# for pos " + str(A_S_L[idx]))
#    plt.plot(X_root[::-1], d_fuel_weight, label = "fuel weight")# for pos " + str(A_S_L[idx]))
#    plt.plot(X_root[::-1], d_strut, label = "strut")# for pos " + str(A_S_L[idx]))
#    plt.plot(X_root[::-1], d_engine, label = "engine")# for pos " + str(A_S_L[idx]))
#    plt.figure()
#    plt.title("Deflection and its derivatives along the half span")

#
    plt.subplot(2, 3, 5)
    plt.plot(X_root, My_dist[idx], label = "My for pos " + str(A_S_L[idx]))
    plt.xlabel("X-position [m]")
    plt.ylabel("My [Nm]")
    plt.title("My distribution")
    ax = plt.gca()
    ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))

#    plt.legend(bbox_to_anchor=(1.05,1), loc="upper left")    
#    plt.subplot(3, 3, 6)
#    plt.plot(label = "strut position" + str(A_S_L[idx]))
#    plt.legend()
#    plt.show()
    
    #
#    plt.subplot(1, 3, 1)
#    plt.plot(X_root[::-1], d, label = "Deflection")# for pos " + str(A_S_L[idx]))
#    plt.title("Deflection [m]")
#    plt.xlabel("X-position [m]")
    
#    plt.subplot(1, 3, 2)
#    plt.plot(X_root[::-1][1:], diff, label = "First derivative")
#    plt.title("First derivative of the deflection")
#    plt.xlabel("X-position [m]")
    
#    plt.subplot(1, 3, 3)
#    plt.plot(X_root[::-1][2:], diff_2, label = "Second derivative")
#    plt.title("Second derivative of the deflection")
#    
#    plt.xlabel("X-position [m]")
#    
#    plt.show()
    
#    
#    

#print(boom_area_all)

