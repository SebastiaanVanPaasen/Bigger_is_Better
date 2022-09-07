import numpy as np
import matplotlib.pyplot as plt
import constants_and_conversions as cc
import centroid_wing as cw
import Airfoil_inertia as ai
import parameter_requirements as pr

from scipy.interpolate import interp1d
from class_I.lift_distr import get_correct_data, lift_distribution, make_avl_file



#D_fus = 7.3

#A_strut = 0.25 * np.pi * ((5 / 1000) ** 2)
#E_strut = 181 * (10 ** 9)  

cr = 5.127
ct = 2.051
b = 9.09
S = 51.518

W_wing = 32965,9463462344
E_wing = 69 * (10 ** 9)
I_zz_wing = 0.0008
L_wing = b / 2

H_cr = 9000
V_cr = 235.42
rho_cr = cc.Rho_0 * ((1 + (cc.a * H_cr) / cc.Temp_0) ** (-(cc.g_0 / (cc.R_gas * cc.a))))
CD_0_cr = 0.0179

W_TO = 1521753.6
W_fuel = 0.
W_N = 0.

#N_eng = 0
#W_eng = 4100 * cc.g_0 + W_N / N_eng
#Rho_fuel = 0.804 * 1000


def deter_d_force(applic, x, force, a_s, I_wing):
    Y_force = (force/ (6 * E_wing * I_wing)) * (2 * ((L_wing - a_s) ** 3) - 3 * ((L_wing - a_s) ** 2) * applic + applic ** 3)
    Theta_force= (-force / (2 * E_wing * I_wing)) * (((L_wing - a_s) - applic) ** 2) * x
    
    Mac = np.array([])
    
    for i in range(len(Theta_force)):
        
        
        if x[i] < applic:
            Mac = np.append(Mac, 0)
        else:
            mac = (force / (6 * E_wing * I_wing[i])) * ((x[i] - applic) ** 3)
            Mac = np.append(Mac, mac)
    

    d_force = Y_force + Theta_force + Mac
    
    return d_force


def get_data(CL):
    make_avl_file()
    
    output_avl = lift_distribution(CL)
    x_pos, cl, cdi = get_correct_data(output_avl)
#    print(x_pos, cl, cdi)

    PolyFitCurveCl = interp1d(x_pos, cl, kind="cubic", fill_value="extrapolate")
    PolyFitCurveidrag = interp1d(x_pos, cdi, kind='cubic', fill_value='extrapolate')
    
    return PolyFitCurveCl,  PolyFitCurveidrag


def calc_chord(x):
    return cr - ((cr - ct) / L_wing) * x
    
    
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
        Volumes = np.append(Volumes, (6.6256 * (10 **(-2))) * (calc_chord(x[i]) ** 2) * width)
        
        
    V_tot = np.sum(Volumes)
    w_spec = w_wing / V_tot


    for i in range(len(x)):
        W = np.append(W, Volumes[i] * w_spec)
        
#    plt.plot(x, W)
#    plt.show()
    
    return W, Volumes


def deter_fuel(w_fuel, volumes, density, x, x_start):
    W = np.array([])
    
    for i in range(len(x)):
        
        
        if x[i] > x_start and w_fuel > 0:
            w_fuel_section = density * volumes[i] * cc.g_0 
        else:
            w_fuel_section = 0


        W = np.append(W, w_fuel_section)
        w_fuel -= w_fuel_section
            
#    plt.plot(x, W)
#    plt.show()
    
    return W
            
    
def indet_sys(F_strut_array, dx, angle, L_s, a_s, a_e, cl_polar, I_wing):   
    X_root = np.arange(0 + dx / 2, L_wing + dx / 2, dx)
    X_tip = np.arange(int(L_wing - a_s) - dx / 2, 0  - dx / 2, -dx)
    x_start = L_wing * 0.34
    
    
    lifts = 3.75 * deter_lift(cl_polar, X_root, dx)  #np.array(len(X_root) * [25295.5]) # 
    drags, thrust = deter_drag(cd_polar, X_root, dx, 3.75)
    weights, Volumes = deter_weight(W_wing, X_root, dx)  #np.array(len(X_root) * [3677]) #
    fuel_weights = deter_fuel(W_fuel / 2, Volumes, Rho_fuel, X_root, x_start)
    

    Lift = []
    Weight = []
    Fuel_weight = []
    
    
#    print("Forces of lift, drag, thrust, weight and fuel weight")
#    print(np.sum(lifts))
#    print(np.sum(drags))
#    print(np.sum(thrust))
#    print(np.sum(weights))
#    print(np.sum(fuel_weights))
#    print()
    

    for i in range(len(X_tip)):
        defl_l = deter_d_force(X_tip[i], np.array([0]), lifts[i], a_s, I_wing)
        defl_w = deter_d_force(X_tip[i], np.array([0]), -weights[i], a_s, I_wing)
        defl_wf = deter_d_force(X_tip[i], np.array([0]), -fuel_weights[i], a_s, I_wing)
           
        Lift = np.append(Lift, defl_l)
        Weight = np.append(Weight, defl_w)
        Fuel_weight = np.append(Fuel_weight, defl_wf)


#    shear_strut = np.sum(lifts[int((L_wing - a_s) / dx):]) - np.sum(weights[int((L_wing - a_s) / dx):]) - np.sum(fuel_weights[int((L_wing - a_s) / dx):])
    mom_lift = np.sum(lifts[int((L_wing - a_s) / dx):] * X_root[:int(a_s / dx)]) - np.sum(weights[int((L_wing - a_s) / dx):] * X_root[:int(a_s / dx)]) - np.sum(fuel_weights[int((L_wing - a_s) / dx):] * X_root[:int(a_s / dx)])
    
#    d_shear = (shear_strut / (6 * E_wing * I_wing)) * (2 * ((L_wing - a_s) ** 3)) 
    d_mom = (mom_lift / (2 * E_wing * I_wing)) * ((L_wing - a_s) ** 2)
    
    d_lift = np.sum(Lift)
    d_weight = np.sum(Weight)
#    d_fuel = np.sum(Fuel_weight)
    
#    print("Deflections of lift, weight and fuel weight front")
#    print(d_lift)
#    print(d_weight)
#    print(d_fuel)
#    print()
#    
#    print("Deflections of shear and moment at the strut")
#    print(d_shear)
#    print(d_mom)
#    print()

        
    d_engine = (-W_eng / (6 * E_wing * I_wing)) * (2 * ((L_wing - a_s) ** 3 ) - 3 * ((L_wing - a_s) ** 2) * (a_e - a_s) + (a_e - a_s) ** 3)
    d_strut_w = (-np.sin(angle) * F_strut_array / (6 * E_wing * I_wing)) * (2 * ((L_wing - a_s) ** 3))
    
#    print("Deflections of engine and strut")
#    print(d_engine)
#    print()
    
    
    d_wing = d_lift + d_weight + d_shear + d_mom + d_engine + d_strut_w + d_fuel
    d_strut = np.sin(angle) * (F_strut_array * L_s) / (E_strut * A_strut)

    diff = abs(d_wing - d_strut)
    idx = np.argmin(diff)
    
#    print("Deflections of wing and strut")
#    print(d_wing[idx])
#    print(d_strut[idx])
#    print()
    
    distributions = [lifts, weights, fuel_weights, W_eng, drags, thrust]
    
    return F_strut_array[idx], d_wing[idx], distributions

    
def strut_opt(A_S_list, A_E, cl_curve, width, I_wing):
    deflections = []
    strut_forces = []
    all_f = []
    L_str = []
    
    for A_S in A_S_list:
        
        L_strut = np.sqrt(D_fus ** 2 + (L_wing - A_S) ** 2)
        gamma = np.arctan(D_fus / (L_wing - A_S))
        print("Angle and length of the strut")
        print(L_strut, gamma)
        
        
        F_strut = np.arange(0, 4000000, 1000)
        force, deflection, all_forces = indet_sys(F_strut, width, gamma, L_strut, A_S, A_E, cl_curve, I_wing)
        
        print("First found optimum")
        print(force, deflection)
        print()
        
        
        F_strut = np.arange(force - 2000, force + 2000, 0.1)
        force, deflection, all_forces = indet_sys(F_strut, width, gamma, L_strut, A_S, A_E, cl_curve, I_wing)
        deflections.append(deflection)
        strut_forces.append(force)
        all_f.append(all_forces)
        L_str.append(L_strut)
        
        print("Final optimum")
        print(force, deflection)
        print()
       
    return strut_forces, deflections, all_f, L_str


A_E = 23
A_S_L = np.arange(5, 21, 5)

cl = W_TO / (0.5 * rho_cr * (V_cr ** 2) * S)
cl_polar, cd_polar = get_data(cl)

dx = 1
results = strut_opt(A_S_L, A_E, cl_polar, dx, I_zz_wing)

X_root = np.arange(0 + dx / 2, L_wing + dx / 2, dx)
X_tip = np.arange(L_wing - dx / 2, 0 - dx / 2, -dx)
X_root_plot = np.append([0], X_root)

Mz_dist = np.zeros((len(A_S_L), len(X_root_plot)))
My_dist = np.zeros((len(A_S_L), len(X_root_plot)))
Vy_dist = np.zeros((len(A_S_L), len(X_root_plot)))
Vz_dist = np.zeros((len(A_S_L), len(X_root_plot)))


for idx in range(len(A_S_L)):
    gamma = np.arctan(D_fus / (L_wing - A_S_L[idx]))
    
    F_str = results[0][idx]
    Lift, Weight, Fuel_weight, W_eng, Drag, Thrust = results[2][idx]

    Lift_mom = Lift * X_root
    Weight_mom = Weight * X_root
    Fuel_mom = Fuel_weight * X_root
#    Eng_mom = W_eng * (L_wing - A_E)
#    Strut_mom = F_str * (L_wing - A_S_L[idx])
    Drag_mom = Drag * X_root
#    Thrust_mom = Thrust * (L_wing - A_E)
    
    Mz_root = sum(Lift_mom) - sum(Weight_mom) - sum(Fuel_mom)# - Eng_mom - Strut_mom
    My_root = sum(Drag_mom)# - Thrust_mom
    Vy_root = sum(Lift) - sum(Weight) - sum(Fuel_weight)# - W_eng - F_str
    Vz_root = Thrust - sum(Drag)
    
    #print("Root moment around z: ", Mz_root)
    print("Root moment around y: ", My_root)
    #print("Root force around y: ", Vy_root)
    print("Root force around z: ", Vz_root)

    
    Mz_dist[idx][0] = Mz_root
    My_dist[idx][0] = My_root
    Vy_dist[idx][0] = Vy_root
    Vz_dist[idx][0] = Vz_root
    
    for i in range(len(X_root_plot) - 1):
        Vy_section = Vy_dist[idx][i] - Lift[i] + Weight[i] + Fuel_weight[i]
        Vz_section = Vz_dist[idx][i] + Drag[i]
        
        if X_root_plot[i] > (L_wing - A_E) and X_root_plot[i - 1] < (L_wing - A_E):
            Vy_section += W_eng
            
        if X_root_plot[i] > (L_wing - A_S_L[idx]) and X_root_plot[i - 1] < (L_wing - A_S_L[idx]):
            Vy_section += F_str
            Vz_section -= Thrust
        
        Vy_dist[idx][i + 1] = Vy_section
        Vz_dist[idx][i + 1] = Vz_section
    
    for i in range(len(X_root_plot) - 1):
        Mz_dist[idx][i + 1] = Mz_dist[idx][i] - Vy_dist[idx][i] * (X_root_plot[i + 1] - X_root_plot[i])
        My_dist[idx][i + 1] = My_dist[idx][i] + Vz_dist[idx][i] * (X_root_plot[i + 1] - X_root_plot[i])
        
    d_lift = 0
    d_weight = 0
    d_fuel_weight = 0
    
    for i in range(len(X_root)):
        d_lift += deter_d_force(X_tip[i], X_root, Lift[i], 0, np.array(len(X_root) * [I_zz_wing])   )
        d_weight += deter_d_force(X_tip[i], X_root, -Weight[i], 0, np.array(len(X_root) * [I_zz_wing]))
        d_fuel_weight += deter_d_force(X_tip[i], X_root, -Fuel_weight[i], 0, np.array(len(X_root) * [I_zz_wing]))
        
    
    d_strut = deter_d_force(A_S_L[idx], X_root, -np.sin(gamma) * F_str, 0, np.array(len(X_root) * [I_zz_wing]))
    d_engine = deter_d_force(A_E, X_root, -W_eng, 0, np.array(len(X_root) * [I_zz_wing]))
    
    
    d = d_lift + d_weight + d_fuel_weight + d_strut + d_engine
    
        
#    plt.subplot(2, 3, 1)
#    plt.plot(X_root_plot, Vy_dist[idx], label = "Vy for pos " + str(A_S_L[idx]))
#    plt.xlabel("X-position [m]")
#    plt.ylabel("Vy [N]")
#    plt.title("Vy distribution")
#    plt.legend()
#    
#    plt.subplot(2, 3, 2)
#    plt.plot(X_root_plot, Vz_dist[idx], label = "Vz for pos " + str(A_S_L[idx]))
#    plt.xlabel("X-position [m]")
#    plt.ylabel("Vz [N]")
#    plt.title("Vz distribution")
#    plt.legend()
#    
#    plt.subplot(2, 3, 3)
##    plt.plot(X_tip, d_lift, label = "lift")
##    plt.plot(X_tip, d_weight, label = "weight")
##    plt.plot(X_tip, d_fuel_weight, label = "fuel weight")
##    plt.plot(X_tip, d_strut, label = "strut")
##    plt.plot(X_tip, d_engine, label = "engine")
#    plt.plot(X_tip, d, label = "Deflection for pos " + str(A_S_L[idx]))
#    plt.xlabel("X-position [m]")
#    plt.ylabel("Deflection [m]")
#    plt.title("Deflection along the span")
#    plt.legend()
#    
#    plt.subplot(2, 3, 4)
#    plt.plot(X_root_plot, Mz_dist[idx], label = "Mz for pos " + str(A_S_L[idx]))
#    plt.xlabel("X-position [m]")
#    plt.ylabel("Mz [Nm]")
#    plt.title("Mz distribution")
#    plt.legend()
#    
#    plt.subplot(2, 3, 5)
#    plt.plot(X_root_plot, My_dist[idx], label = "My for pos " + str(A_S_L[idx]))
#    plt.xlabel("X-position [m]")
#    plt.ylabel("My [Nm]")
#    plt.title("My distribution")
#    plt.legend()
#    
#    plt.show()



l_spar_h, t_spar_v, t_spar_h = cw.l_spar_h, cw.t_spar_v, cw.t_spar_h
N = L_wing / dx

Izz_list = np.zeros((len(A_S_L), len(X_root_plot)))
boom_area = np.zeros((len(A_S_L), len(X_root_plot)))
Mz_total = np.zeros((len(A_S_L), len(X_root_plot)))


for idx in range(len(A_S_L)):
    I_zz_spar, I_yy_spar, I_yz_spar = ai.I_zz_spars(l_spar_h, t_spar_v, t_spar_h, N, b ,calc_chord)
    I_zz_req = pr.required_Izz(N, b, calc_chord, Mz_dist[idx][1:])
#    print(I_zz_req)
    airfoil_area, z_c_airfoil, y_c_airfoil = cw.get_skin_centroid(b, N, calc_chord)
    boom_area = ai.wing_geometry(I_zz_req, I_zz_spar, N, b, calc_chord)[0][0]
#    print(boom_area)
    I_zz_wing, I_yy_wing, I_yz_wing = ai.inertia_wing(I_zz_spar, I_yy_spar, I_yz_spar, boom_area, N, b, calc_chord)

    results = strut_opt(A_S_L, A_E, cl_polar, dx, I_zz_wing[int((L_wing - A_S_L[idx]) / dx)])
    
    gamma = np.arctan(D_fus / (L_wing - A_S_L[idx]))
    

    F_str = results[0][idx]
    Lift, Weight, Fuel_weight, W_eng, Drag, Thrust = results[2][idx]
    
    
    Lift_mom = Lift * X_root
    Weight_mom = Weight * X_root
    Fuel_mom = Fuel_weight * X_root
    Eng_mom = W_eng * (L_wing - A_E)
    Strut_mom = F_str * (L_wing - A_S_L[idx])
    
    Mz_root = sum(Lift_mom) - sum(Weight_mom) - sum(Fuel_mom) - Eng_mom - Strut_mom
    Vy_root = sum(Lift) - sum(Weight) - sum(Fuel_weight) - W_eng - F_str
    
#    print(Mz_root)
#    print(Vy_root)
    
    
    Mz_dist[idx][0] = Mz_root
    Vy_dist[idx][0] = Vy_root
#    D[idx][0] = 0
    
    for i in range(len(X_root_plot) - 1):
        Vy_section = Vy_dist[idx][i] - Lift[i] + Weight[i] + Fuel_weight[i]
        
        if X_root_plot[i] > (L_wing - A_E) and X_root_plot[i - 1] < (L_wing - A_E):
            Vy_section += W_eng
            
        if X_root_plot[i] > (L_wing - A_S_L[idx]) and X_root_plot[i - 1] < (L_wing - A_S_L[idx]):
            Vy_section += F_str
            
        Vy_dist[idx][i + 1] = Vy_section
        
    
    for i in range(len(X_root_plot) - 1):
        Mz_dist[idx][i + 1] = Mz_dist[idx][i] - Vy_dist[idx][i] * (X_root_plot[i + 1] - X_root_plot[i])
#        D[idx][i + 1] = -(Mz_dist[idx][i + 1] * (X_root_plot[i + 1]) ** 2) / (2 * E_wing * I_wing)
    
    
    d_lift = 0
    d_weight = 0
    d_fuel_weight = 0
    for i in range(len(X_root)):
    
        d_lift += deter_d_force(X_tip[i], X_root, Lift[i], 0, I_zz_wing)   
        d_weight += deter_d_force(X_tip[i], X_root, -Weight[i], 0, I_zz_wing)
        d_fuel_weight += deter_d_force(X_tip[i], X_root, -Fuel_weight[i], 0, I_zz_wing)
        
        
#    d_strut = deter_d_force(A_S_L[idx], X_root, -np.sin(gamma) * F_str, 0, I_zz_wing)
#    d_engine = deter_d_force(A_E, X_root, -W_eng, 0, I_zz_wing)
    
    
#    d = d_lift + d_weight + d_fuel_weight + d_strut + d_engine

L_strut_list = results[3]
print(L_strut_list)
#    print(d)
#    print(len(d))
#    plt.subplot(1, 3, 1)
#    plt.plot(X_root_plot, Vy_dist[idx], label = "Vy for pos " + str(A_S_L[idx]))
#    plt.xlabel("X-position 1 [m]")
#    plt.ylabel("Vy [N]")
#    plt.title("Fy distribution")
#    plt.legend()
#    
#    plt.subplot(1, 3, 2)
#    plt.plot(X_root_plot, Mz_dist[idx], label = "Mz for pos " + str(A_S_L[idx]))
#    plt.xlabel("X-position 2 [m]")
#    plt.ylabel("Mz [Nm]")
#    plt.title("Mz distribution")
#    plt.legend()
#
#    plt.subplot(1, 3, 3)
##    plt.plot(X_tip, d_lift, label = "lift")
##    plt.plot(X_tip, d_weight, label = "weight")
##    plt.plot(X_tip, d_fuel_weight, label = "fuel weight")
##    plt.plot(X_tip, d_strut, label = "strut")
##    plt.plot(X_tip, d_engine, label = "engine")
#    plt.plot(X_tip, d[idx], label = "Deflection for pos " + str(A_S_L[idx]))
#    plt.xlabel("X-position 3 [m]")
#    plt.ylabel("Deflection [m]")
#    plt.title("Deflection along the span")
#    plt.legend()
#    
##    plt.subplot(2, 2, 4)
##    plt.plot(X_root_plot, D[idx], label = "D based on Mz for pos " + str(A_S_L[idx]))
#    
#    plt.show()



