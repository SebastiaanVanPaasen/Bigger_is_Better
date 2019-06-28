import numpy as np
import matplotlib.pyplot as plt
import constants_and_conversions as cc
import centroid_wing as cw
import Airfoil_inertia as ai
import parameter_requirements as pr
import math as m

from scipy.interpolate import interp1d
from class_I.lift_distr import get_correct_data, lift_distribution

D_fus = 7.3

#R_strut = 5 / 1000
#A_strut = 0.25 * np.pi * ((2 * R_strut) ** 2) ##
E_strut = 130 * (10 ** 9)  
AR = 13
taper = 0.357 
MAC = 4.3#4.03
sigma_strut = 1400 * (10 ** 6)
density_strut = 1580.

cr = 5.89
ct = 2.1
b = 52
S = 208
dihedral = (1.5/180) * np.pi

LE_root = 19.11
pos_of_chord = 0.25
sweep = 0
pos_from_centerline = 2.5
strut_loc_fus = LE_root + pos_of_chord * cr + np.tan(sweep)* (b/2)
strut_heigth = D_fus - np.tan(dihedral)*(b/2) - 1.

W_wing = 25480*9.81/2 #139716/2
E_wing = 75 * (10 ** 9)
I_zz_wing = 0.25
L_wing = b / 2

H_cr = 9000#7636.9
V_cr = 201#272.5
rho_cr = cc.Rho_0 * ((1 + (cc.a * H_cr) / cc.Temp_0) ** (-(cc.g_0 / (cc.R_gas * cc.a) + 1)))
CD_0_cr = 0.021323
#print(rho_cr)
W_TO = 1497151
W_fuel = 170548
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
        Volumes = np.append(Volumes, (0.0809 * (calc_chord(x[i]) ** 2) * width))
        
        
    V_tot = np.sum(Volumes)
    w_spec = w_wing / V_tot
    

    for i in range(len(x)):
        W = np.append(W, Volumes[i] * w_spec)
        
#    plt.plot(x, W)
#    plt.show()
    
    return W, Volumes


def deter_fuel(w_fuel, volumes, density, x, x_start):
    W = np.array([])
#    print("total volume", sum(volumes))
#    print(w_fuel)
    for i in range(len(x)):
        
        
        if x[i] > x_start and w_fuel > 0:
            w_fuel_section = density * volumes[i]* 0.6 * cc.g_0 * 0.98
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
    
    alpha = np.radians(2.235)
    n = 2.56 * 1.5
    
    
    lifts_v = n * deter_lift(cl_polar, X_tip[::-1], dx)   
#    print("weight", W_TO)
#    print("total lift", 2*sum(lifts_v))
    drags_v, thrust = deter_drag(cd_polar, X_tip[::-1], dx, n)   
    weights_v, Volumes = deter_weight(W_wing, X_tip[::-1], dx)  #
    fuel_weights_v = deter_fuel(W_fuel / 2, Volumes, Rho_fuel, X_tip[::-1], x_start)
    
    lifts = np.cos(alpha) * lifts_v + np.sin(alpha) * drags_v
    drags = np.sin(alpha) * (-lifts_v + weights_v + fuel_weights_v + W_eng_v) + np.cos(alpha) * drags_v 
    
    weights = np.cos(alpha) *  weights_v
    fuel_weights = np.cos(alpha) * fuel_weights_v
    W_eng = np.cos(alpha) * W_eng_v
    

#    print("Forces of lift, drag, thrust, weight and fuel weight")
#    print(np.sum(np.sin(alpha) * (-lifts_v)))
#    print(np.sum(np.sin(alpha) * ( weights_v + fuel_weights_v + W_eng_v)))
#    print(np.sum(np.cos(alpha) * drags_v))
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
    d_strut = np.sin(angle) * sigma_strut * L_s / E_strut 
#    print("deflection of strut", d_strut)
    d_strut_v = d_wing[int((a_s) / dx)] - d_strut
#    print("dstrutv",d_strut_v)
    y_a = -(2 * L_wing ** 3 - 3 * (L_wing ** 2) *  (a_s + (dx/2)) +  (a_s + (dx/2)) ** 3)/(6 * E_wing * I_zz_sections[::-1][int((a_s)/dx)])
    theta_a = (((L_wing - (a_s + dx/2)) ** 2) / (2 * E_wing * I_zz_sections[::-1][int((a_s)/dx)])) * (a_s+(dx/2))
#    print("y_a", y_a)
#    print("theta", theta_a)
    
    F_strut = -(d_strut_v / (y_a + theta_a))  # (2 * ((L_wing) ** 3))) * (6 * E_wing * I_wing)) / np.sin(angle)
    
    distributions = [lifts, weights, fuel_weights, W_eng, drags, thrust]
    
    return F_strut, d_strut, d_strut_v, d_wing, distributions



A_E = 19.65 #7m from center of fuselage
A_S_L = np.arange(8.05, 8.55, 0.5)

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
#boom_area_old = 0.5
boom_area_all = np.zeros(len(A_S_L)) 
 
F_strut = np.zeros(len(A_S_L))
L_str = np.zeros(len(A_S_L))
gamma_all = np.zeros(len(A_S_L))
R_strut_new = np.zeros(len(A_S_L))
A_req_new = np.zeros(len(A_S_L))
strut_mass = np.zeros(len(A_S_L))

for idx in range(len(A_S_L)):
    boom_area_old = np.ones(len(X_root))*0.5
    I_zz_spar, I_yy_spar, I_yz_spar = ai.I_zz_spars(l_spar_h, t_spar_v, t_spar_h, N, b ,calc_chord, boom_area_old,X_root, dx)
    I_zz_req = pr.required_Izz(N, b, calc_chord, Mz_dist[idx], boom_area_old, X_root, dx)

    boom_area_new = ai.wing_geometry(I_zz_req, I_zz_spar, N, b, calc_chord, boom_area_old, X_root, dx)

#    print(abs(boom_area_new - boom_area_old) )
    z = 0
    while abs(boom_area_new[0] - boom_area_old[0]) > 1 / 100000 or z < 5:
#        print(abs(boom_area_new - boom_area_old) )
        print("New iteration", z)
        boom_area_old = boom_area_new
#        boom_area_all = 0.0045
        
        I_zz_spar, I_yy_spar, I_yz_spar = ai.I_zz_spars(l_spar_h, t_spar_v, t_spar_h, N, b ,calc_chord, boom_area_old, X_root, dx)
        I_zz_req = pr.required_Izz(N, b, calc_chord, Mz_dist[idx], boom_area_old, X_root, dx)
        
        boom_area_new = ai.wing_geometry(I_zz_req, I_zz_spar, N, b, calc_chord, boom_area_old,X_root, dx)
        airfoil_area, z_c_airfoil, y_c_airfoil = cw.get_skin_centroid(N, b, calc_chord,dx, cw.t_skin)
        
        z_centroid_all_sec, y_centroid_all_sec, y_loc_spar_up, y_loc_spar_low, y_loc_stiff_up, y_loc_stiff_low, y_vertical_spar, z_loc_stiff_up, spar_loc_sec, z_loc_stiff_low, spar_areas_verti = cw.wing_centroid(boom_area_new, cw.spar_areas_hori, cw.t_spar_v, z_c_airfoil, y_c_airfoil, cw.n_stiff_up, cw.n_stiff_low, N, b, calc_chord, X_root, dx)
        #
        #    
        #    #    print("I_zz_req",I_zz_req)
#        boom_area = boom_area_new
        boom_area_all[idx] = max(boom_area_new) * 10000
        
        print("Updated boom area for strut pos" + str(A_S_L[idx]))
#        
        print(max(boom_area_new) * 10000)
#        print()
##        
    
        I_zz_sections, I_yy_wing, I_yz_wing = ai.inertia_wing(I_zz_spar, I_yy_spar, I_yz_spar, boom_area_new, N, b, calc_chord, X_root, dx)
    #       
#        print("I_zz_sections",I_zz_sections)
#        print(I_yy_wing)
#        print(I_yz_wing)
    #        I_zz_sections = np.array(293 * [0.00499])
        strut_loc_wing = L_wing - A_S_L[idx] - pos_from_centerline
           
        gamma = np.arctan(strut_heigth / strut_loc_wing)
        gamma_all[idx]=gamma
        L_strut = np.sqrt(strut_heigth**2 + strut_loc_wing**2)
        L_str[idx] = L_strut
        print("length strut", L_strut)
        
        
        F_str, d_str, d_str_v, d_w, all_forces = indet_sys(dx, gamma, L_strut, A_S_L[idx], A_E, cl_polar, I_zz_sections)
        print("Strut force and deflections")
        print("Strut force", F_str)
        print("Deflection of the strut", d_str)
        print(d_str_v)
#        print("Required strut area", F_str / (140 * (10 ** 6)))
        #        print()
        F_strut[idx] = F_str
        F_str = 0.
        #    F_str = results[0]
        #    print("Strut force in first optimisation")
        #    print(F_str)
        #    print()
        
        Lift, Weight, Fuel_weight, W_eng, Drag, Thrust = all_forces
        
        Lift_mom = Lift * X_tip[::-1]
        Weight_mom = Weight * X_tip[::-1]
        Fuel_mom = Fuel_weight * X_tip[::-1]
        Eng_mom = W_eng * (L_wing - A_E)
        Strut_mom = F_str * (L_wing - A_S_L[idx])
        #        Strut_mom_hori = F_str/np.tan(angle) * y_centroid_all_sec[int((L_wing - A_S_L[idx])/dx)] - 
        Drag_mom = Drag * X_tip[::-1]
        Thrust_mom = Thrust * (L_wing - A_E)
        
        Strut_z = ((F_str/np.tan(gamma)) * (-min(y_loc_stiff_low[int((L_wing - A_S_L[idx])/dx)]) + y_centroid_all_sec[0]))
        print(Strut_z)
        Strut_y = ((F_str/np.tan(gamma)) * (-0.25*calc_chord(X_root[0]) + z_centroid_all_sec[0]))
        print(Strut_y)
        
        Mz_root = sum(Lift_mom) - sum(Weight_mom) - sum(Fuel_mom) - Eng_mom - Strut_mom - Strut_z
        My_root = sum(Drag_mom) - Thrust_mom - Strut_y
        Vy_root = sum(Lift) - sum(Weight) - sum(Fuel_weight) - W_eng - F_str
        Vz_root = Thrust - sum(Drag)
        Vx_root = (F_str/np.tan(gamma))
        
        print("Root moment around z: ", Mz_root)
        print("Root moment around y: ", My_root)
        print("Root force around y: ", Vy_root)
        print("Root force around z: ", Vz_root)
        print("Root force around x: ", Vx_root)
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
            Mz_dist[idx][i + 1] = Mz_root - Vy_root * X_root[i+1] + Vx_root * (-y_centroid_all_sec[int(X_root[i+1]/dx)] + y_centroid_all_sec[0])  + sum(((Lift - Weight - Fuel_weight) * (X_root[i+1] - X_tip[::-1]))[:i])
            My_dist[idx][i + 1] = My_root + Vz_root * X_root[i+1] + Vx_root * (z_centroid_all_sec[0] - (z_centroid_all_sec[int(X_root[i + 1]/dx)] - 0.25 * calc_chord(X_root[i + 1]) + 0.25 * calc_chord(0))) + sum((Drag*(X_root[i+1] - X_tip[::-1]))[:i])
            
            if X_root[i+1] > (L_wing - A_E): #and X_root[i - 1] < (L_wing - A_E):
                Mz_dist[idx][i + 1] -= W_eng * (X_root[i+1]- (L_wing - A_E))
                
            if X_root[i+1] > (L_wing - A_S_L[idx]): #and X_root[i - 1] < (L_wing - A_S_L[idx]):
                Mz_dist[idx][i + 1] -= F_str * (X_root[i+1]- (L_wing - A_S_L[idx]))
                Mz_dist[idx][i + 1] += (F_str/np.tan(gamma)) * (-min(y_loc_stiff_low[int((L_wing - A_S_L[idx])/dx)]) + y_centroid_all_sec[int(X_root[i+1]/dx)])#((-min(y_loc_stiff_low[int((L_wing - A_S_L[idx])/dx)]) - y_centroid_all_sec[int(X_root[i + 1]/dx)]) - (-min(y_loc_stiff_low[int((L_wing - A_S_L[idx])/dx)]) - y_centroid_all_sec[int((L_wing - A_S_L[idx])/dx)]))
                My_dist[idx][i + 1] += (F_str/np.tan(gamma)) * (z_centroid_all_sec[int(X_root[i + 1]/dx)] - 0.25 * calc_chord(X_root[i + 1]))# + (0.25 * calc_chord(X_root[int((L_wing - A_S_L[idx])/dx)]) - 0.25 * calc_chord(X_root[i + 1])))#(z_centroid_all_sec[int(X_root[i+1]/dx)]) + (-0.25*calc_chord(L_wing - A_S_L[idx]) + z_centroid_all_sec[int((L_wing - A_S_L[idx])/dx)]))
#                print(" distance" )
#                print((-(-0.25*calc_chord(L_wing - A_S_L[idx]) + z_centroid_all_sec[int(X_root[i+1]/dx)]), (-0.25*calc_chord(L_wing - A_S_L[idx]) + z_centroid_all_sec[int((L_wing - A_S_L[idx])/dx)])))
            
            if X_root[i+1] > (L_wing - A_E): #and X_root[i - 1] < (L_wing - A_E):
        
                My_dist[idx][i + 1] -= Thrust * (X_root[i+1]- (L_wing - A_E))
        z+=1
    
#    boom_area_all.append(boom_area_new)
    
    K = 1
    
    I_req_buck = (abs(F_strut[idx])/np.sin(gamma_all[idx]))*(K*L_str[idx])**2/(np.pi**2*E_strut)
    R_strut_new[idx] = (I_req_buck/(0.25*np.pi))**0.25 
    A_req_new[idx] = abs(F_strut[idx])/sigma_strut

#    print("r_strut_buck", R_strut_new)
#    print("r_strut_tens", np.sqrt(A_req_new/np.pi))
    
    strut_volume = np.pi*(R_strut_new[idx])**2*L_str[idx]
    strut_mass[idx] = strut_volume*density_strut
    
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
    
#    d_strut = 0

    d = d_lift + d_weight + d_fuel_weight + d_strut + d_engine

#    x_start = 0*(b/2)
#    plt.figure()
#    plt.plot(X_tip[::-1], all_forces[2], label = "start fuel tank " + str(round(x_start,2)))
#    plt.xlabel("X-position [m]")
#    plt.ylabel("Fuel weigth [N]")
#    plt.show()
#    plt.legend()
#        
#    for i in range(1, len(d)):
#        diff.append(d[i - 1] - d[i])
#        
# 
#    for i in range(1, len(diff)):
#        diff_2.append(diff[i - 1] - diff[i])
            

#    print("deflection of components after strut calculation")
#    print(d_lift[int((A_S_L[idx])/dx)])
#    print(d_weight[int((A_S_L[idx])/dx)])
#    print(d_fuel_weight[int((A_S_L[idx])/dx)])
#    print(d_engine[int((A_S_L[idx])/dx)])
#    print(d_strut[int((A_S_L[idx])/dx)])
#    print(d[int((A_S_L[idx])/dx)])
#    print("izz", I_zz_sections[::-1])
#    print(Vy_dist[idx][0])
#    print(Mz_dist[idx][0])
#    print(My_dist[idx][0])
#    print(Vz_dist[idx][0])
#    plt.rcParams.update({'font.size': 15})        
    plt.figure()
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
    
    label_list = [21, 19, 17, 15, 13, 11, 9]
    plt.subplot(2, 3, 5)
    plt.plot(X_root, My_dist[idx], label = "Strut position " + str(label_list[idx]))
    plt.xlabel("X-position [m]")
    plt.ylabel("My [Nm]")
    plt.title("My distribution")
    ax = plt.gca()
    ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
##
#    plt.legend(bbox_to_anchor=(1.05,1), loc="upper left")    
#    plt.subplot(3, 3, 6)
#    plt.plot(label = "strut position" + str(A_S_L[idx]))
#    plt.legend()
#    plt.show()
    
print("r_strut_buck", R_strut_new)
print("r_strut_tens", np.sqrt(A_req_new/np.pi))
print("mass_strut", strut_mass)

#    plt.subplot(1, 3, 1)
#    plt.plot(X_root[::-1], d, label = "Deflection")# for pos " + str(A_S_L[idx]))
#    plt.title("Deflection [m]")
#    plt.xlabel("X-position [m]")
#    
#    plt.subplot(1, 3, 2)
#    plt.plot(X_root[::-1][1:], diff, label = "First derivative")
#    plt.title("First derivative of the deflection")
#    plt.xlabel("X-position [m]")
#    
#    plt.subplot(1, 3, 3)
#    plt.plot(X_root[::-1][2:], diff_2, label = "Second derivative")
#    plt.title("Second derivative of the deflection")
#    plt.xlabel("X-position [m]")
##    
#    plt.show()
    
#def strut_area_req_F(sigma):
#    
#    F_strut = prac.F_strut
#    A_req = np.zeros(len(prac.A_S_L))
#    
#    for i in range(len(prac.A_S_L)):
#
#        A_req[i] = abs(F_strut[i])/sigma
#
#    return A_req
##
#def strut_area_req_B(F_strut, L_strut, gamma, A_S_L, E_strut):
#    
#    P_cr = abs(F_strut)
#    I_req = np.zeros(len(A_S_L))
#    req_area = np.zeros(len(A_S_L))
#    sigma_crit = np.zeros(len(A_S_L))
#    K = 1. 
#    
#    for i in range(len(A_S_L)):
#        I_req[i] = (P_cr[i]/np.sin(gamma[i]))*(K*L_strut[i])**2/(np.pi**2*E_strut)
#
#        R_strut_new = (I_req[i]/(0.25*np.pi))**0.25
##        while I_strut <I_req[i]:
##            R_strut += 0.001
##            I_strut = 0.25*np.pi*(R_strut**4)
##        print(I_req[i])
#
##        print("r_strut", R_strut_new)
#
#        req_area[i] = np.pi*R_strut_new**2
#        sigma_crit[i] = P_cr[i]/req_area[i]
#        
#    return I_req, req_area, sigma_crit
#
#def strut_cost(A_req_B, density, cost, L_strut, A_S_L):
#    Max_area = np.zeros(len(A_S_L))
#    strut_volume = np.zeros(len(A_S_L))
#    strut_mass = np.zeros(len(A_S_L))
#    strut_cost = np.zeros(len(A_S_L))
#
#    
#    for i in range(len(A_S_L)):
#        Max_area[i] = A_req_B[i]#max(A_req_P[i], A_req_B[i])#A_req_P[i]#
#        
#        strut_volume[i] = Max_area[i]*L_strut[i]
#        strut_mass[i] = strut_volume[i]*density
#        strut_cost[i] = strut_mass[i]*cost
#
#    return Max_area, strut_volume, strut_mass, strut_cost
#
#
#def wing_price_weight(A_req_B, density, cost, N, t_skin, b, qcsweep, dx, L_strut, A_S_L, boom_area):
#            
#    Max_area, strut_volume, strut_mass, cost_strut = strut_cost(A_req_B, density, cost, L_strut, A_S_L)
#    airfoil_area, z_c_airfoil, y_c_airfoil = cw.get_skin_centroid(N, b, calc_chord, dx, t_skin)
##    boom_area = boom_area_all
##    print("boom area",boom_area)
#    X_root = np.arange(0, (b/2)+dx, dx)
#
#    spar_areas_verti = cw.wing_centroid(boom_area, cw.spar_areas_hori, cw.t_spar_v, z_c_airfoil, y_c_airfoil, cw.n_stiff_up, cw.n_stiff_low, N, b, calc_chord,X_root, dx)[10]
#    spar_areas_hori = cw.spar_areas_hori
##    print(spar_areas_hori)
#    nr_stiff  = cw.n_stiff_low + cw.n_stiff_up
##    Sweep_LE = m.atan(m.tan(qcsweep) - 4 / AR * (-0.25 * (1 - taper) / (1 + taper))) # rad
#    
##    print(spar_areas_verti[0])
##    print(spar_areas_hori)
##    print(spar_areas_hori)
##    print(spar_areas_verti[0])
##    print(spar_areas_verti[1])
##    print(len(spar_areas_verti[0]))
#    
#    total_spar_volume = 0
#    
#    for i in range(len(spar_areas_verti)):
#        spar_volume = 0
#        for j in range(len(spar_areas_verti[0])):
#            
#            spar_volume += (spar_areas_verti[i][j] + spar_areas_hori[i]*2) 
#        
#        total_spar_volume += spar_volume*dx
#        
#    total_boom_volume = np.zeros(len(A_S_L))
#
#    for i in range(len(A_S_L)):
#        
#        for j in range(len(X_root)):
#            
#            if boom_area[j]> 0:
#            
#                total_boom_volume[i] += (boom_area[j] * nr_stiff)*dx
#        
#            else:
#                total_boom_volume[i] += 0
#        
##    print("boom_area", boom_area)
#    boom_mass = total_boom_volume * density
#    boom_cost = boom_mass * cost
#    
#    skin_volume = np.zeros(len(X_root))
#
#    for i in range(len(X_root)-1):
#        
#        skin_volume[i] = ai.s_airfoil(N,b,calc_chord, X_root)[0][i] * dx * t_skin[i]
##    print(skin_volume)
#    skin_mass = sum(skin_volume) * density
#    skin_price = skin_mass * cost
#    
#    spar_mass = total_spar_volume * density
#    spar_price = spar_mass * cost
#    
#    total_price = np.zeros(len(A_S_L))
#    total_mass = np.zeros(len(A_S_L))
##    print(spar_length, spar_volume, spar_mass, total_spar_area)
##    print(skin_price, spar_price)
#    for i in range(len(A_S_L)):
#        
#        total_price[i] = (skin_price + spar_price)*2 + cost_strut[i]*2 + boom_cost[i]*2
#        total_mass[i] = (skin_mass + spar_mass)*2 + strut_mass[i]*2 + boom_mass[i]*2
#    
#    return  skin_mass, spar_mass, boom_mass, total_mass, total_price, boom_cost
##
#
#
#
#
#
#t_skin = cw.t_skin
#density_wing = 2780
#cost_wing = 1.96 
#density_strut = 1580
#cost_strut = 35.2
#qcsweep = 0 * np.pi
#
#
#I_req, A_req_B, sigma_crit= strut_area_req_B(F_strut, L_str, gamma_all, A_S_L, E_strut)
#
#max_strut_area, strut_volume, strut_mass, cost_strut = strut_cost(A_req_B, density_strut, cost_strut, L_str, A_S_L)
#skin_mass, spar_mass, boom_mass, total_mass, total_price, boom_cost = wing_price_weight(A_req_B, density_wing, cost_wing, N, t_skin, b, qcsweep, dx, L_str, A_S_L, boom_area_new)
#
#
#print("area due to buckling", A_req_B)
#print()
#print("Total spar mass", spar_mass*2)
#print()
##print("critical stress", sigma_crit)
##print()
#print("Total skin mass", skin_mass*2)
#print()
##print("spar mass", spar_mass)
##print()
#print("Total boom mass", boom_mass*2)
#print()
#print("Total boom cost", boom_cost*2)
#print()
#print("Total strut mass", strut_mass*2)
#print()
#print("Total strut cost", cost_strut*2)
#print()
#print("Total wing mass", total_mass)
#print()
#print("Total wing cost", total_price)
#
