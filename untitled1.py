import numpy as np
import matplotlib.pyplot as plt
import constants_and_conversions as cc

from scipy.interpolate import interp1d
from class_I.lift_distr import get_correct_data, lift_distribution, make_avl_file

D_fus = 7.3

A_strut = 0.25 * np.pi * ((5 / 1000) ** 2)
E_strut = 181 * (10 ** 9)  

cr = 6.04
ct = 1.89
b = 60
S = 235.7

W_wing = 200737.6
E_wing = 69 * (10 ** 9)
I_wing = 0.000499
L_wing = b / 2

H_cr = 9000
V_cr = 2335.42
rho_cr = cc.Rho_0 * ((1 + (cc.a * H_cr) / cc.Temp_0) ** (-(cc.g_0 / (cc.R_gas * cc.a))))

W_TO = 1521753.6
W_fuel = 200716.1
W_eng = 73144  #4100 * cc.g_0
Rho_fuel = 0.804 * 1000


def deter_d_distr(applic, x, distr):
    Y_distr = (distr / (24 * E_wing * I_wing)) * ((L_wing - applic) ** 3) * (3 * L_wing + applic)
    Theta_distr = (-distr / (6 * E_wing * I_wing)) * ((L_wing - applic) ** 3) * x
    Mac = np.array([])
    for i in range(len(Theta_distr)):
        if x[i] < applic:
            Mac = np.append(Mac, 0)
        else:
            mac = (distr / (24 * E_wing * I_wing)) * ((x[i] - applic) ** 4)
            Mac = np.append(Mac, mac)
        
    d_distr= Y_distr+ Theta_distr + Mac
    
    return d_distr


def deter_d_force(applic, x, force, a_s):
    Y_force = (force/ (6 * E_wing * I_wing)) * (2 * ((L_wing - a_s) ** 3) - 3 * ((L_wing - a_s) ** 2) * applic + applic ** 3)
    Theta_force= (-force / (2 * E_wing * I_wing)) * (((L_wing - a_s) - applic) ** 2) * x
    
    Mac = np.array([])
    
    for i in range(len(Theta_force)):
        
        
        if x[i] < applic:
            Mac = np.append(Mac, 0)
        else:
            mac = (force / (6 * E_wing * I_wing)) * ((x[i] - applic) ** 3)
            Mac = np.append(Mac, mac)
    

    d_force = Y_force + Theta_force + Mac
    
    return d_force

def deter_strut(F_strut_array, angle, L_s, applications, forces):
    a_l, a_w, a_eng, a_s = applications
    L_dist, W_dist, W_eng = forces
    
    Lift = L_dist * ((L_wing - a_s) ** 4) / (8 * E_wing * I_wing)
    
    if a_l < a_s:
        Strut_shear = (L_dist * (a_s - a_l)) * ((L_wing - a_s) ** 3) / (3 * E_wing * I_wing)
        Strut_moment = (L_dist * ((a_s - a_l) ** 2) / 2) * ((L_wing - a_s) ** 2) / (2 * E_wing * I_wing)
    else:
        Strut_shear = 0
        Strut_moment = 0
        
        
    Weight = -W_dist * ((L_wing - a_s) ** 4) / (8 * E_wing * I_wing)
    
    if a_w < a_s:
        Strut_shear += (-W_dist * (a_s - a_w)) * ((L_wing - a_s) ** 3) / (3 * E_wing * I_wing)
        Strut_moment += (-W_dist * ((a_s - a_w) ** 2) / 2) * ((L_wing - a_s) ** 2) / (2 * E_wing * I_wing)

    
    Engine = (-W_eng / (6 * E_wing * I_wing)) * (2 * ((L_wing - a_s) ** 3 ) - 3 * ((L_wing - a_s) ** 2) * (a_eng - a_s) + (a_eng - a_s) ** 3)

    
    Strut = (-np.sin(angle) * F_strut_array * ((L_wing - a_s) ** 3)) / (3 * E_wing * I_wing)
    
    
    d_wing = Lift + Weight + Strut_shear + Strut_moment + Engine + Strut 
    d_strut = np.sin(angle) * (F_strut_array * L_s) / (E_strut * A_strut)
    
    diff = abs(d_wing - d_strut)
    idx = np.argmin(diff)
    
#    print(Lift, Weight, Strut_shear, Strut_moment, Engine)
    return F_strut_array[idx]


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
            
    
    

#l_distr = 25295.5
#w_distr = 3677
#
#force = l_distr, w_distr, W_eng

A_L = 0
A_W = 0
A_E = 23
A_S = 5

applic = A_L, A_W, A_E, A_S

L_strut = np.sqrt(D_fus ** 2 + (L_wing - A_S) ** 2)
gamma = np.arctan(D_fus / (L_wing - A_S))
#print(L_strut, gamma)


def indet_sys(F_strut_array, dx, angle, L_s, a_s, a_e):
    cl = W_TO / (0.5 * rho_cr * (V_cr ** 2) * S)
    cl_polar, cd_curve = get_data(cl)
    
    
    X_root = np.arange(0 + dx / 2, int(L_wing - a_s) + dx / 2, dx)
    X_tip = np.arange(int(L_wing - a_s), 0  + dx / 2, -dx)
#    x_start = L_wing * 0.35
    print(int(L_wing - a_s))
    
    lifts = np.array(len(X_root) * [25295.5]) #deter_lift(cl_polar, X, dx) 
    weights = np.array(len(X_root) * [3677])  #, Volumes = deter_weight(W_wing, X, dx)
#    fuel_weights = deter_fuel(W_fuel / 2, Volumes, Rho_fuel, X, x_start)
    
    
    Lift = []
    Weight = []
#    Fuel_weight = []
    

    for i in range(len(X_tip)):
        defl_l = deter_d_force(X_tip[i], np.array([0]), lifts[i], a_s)
        defl_w = deter_d_force(X_tip[i], np.array([0]), -weights[i], a_s)
#        defl_wf = deter_d_force(X_tip[i], np.array([0]), -fuel_weights[i], a_s)
           
        Lift = np.append(Lift, defl_l)
        Weight = np.append(Weight, defl_w)
#        Fuel_weight = np.append(Fuel_weight, defl_wf)

    shear_strut = np.sum(lifts[int(L_wing - a_s - a_s):]) - np.sum(weights[int(L_wing - a_s - a_s):]) #- np.sum(fuel_weights[int(L_wing - a_s - a_s):])
    mom_lift = np.sum(lifts[int(L_wing - a_s - a_s):] * X_root[:a_s]) - np.sum(weights[int(L_wing - a_s - a_s):] * X_root[:a_s]) # np.sum(fuel_weights[int(L_wing - a_s - a_s):] * X_root[:a_s])
    
    d_shear = (shear_strut / (6 * E_wing * I_wing)) * (2 * ((L_wing - a_s) ** 3)) # - 3 * (L_wing ** 2) * a_s + a_s ** 3)
    d_mom = (mom_lift / (2 * E_wing * I_wing)) * ((L_wing - a_s) ** 2)
    
    d_lift = np.sum(Lift)
    d_weight = np.sum(Weight)
#    d_fuel = np.sum(Fuel_weight)
    
    print("Deflections of lift, weight and fuel weight front")
    print(d_lift)
    print(d_weight)
#    print(d_fuel)
    print()
    
    print("Deflections of shear and moment at the strut")
    print(d_shear)
    print(d_mom)
    print()

        
    d_engine = (-W_eng / (6 * E_wing * I_wing)) * (2 * ((L_wing - a_s) ** 3 ) - 3 * ((L_wing - a_s) ** 2) * (a_e - a_s) + (a_e - a_s) ** 3)
    d_strut_w = (-np.sin(angle) * F_strut_array / (6 * E_wing * I_wing)) * (2 * ((L_wing - a_s) ** 3))
    
    print("Deflections of engine and strut")
    print(d_engine)
    print(d_strut_w)
    print()
    
    
    d_wing = d_lift + d_weight + d_shear + d_mom + d_engine + d_strut_w #+ d_fuel
    print(d_wing)
    d_strut = np.sin(angle) * (F_strut_array * L_s) / (E_strut * A_strut)
    print(d_strut)
#    print(d_strut)
    print(d_wing[1086])
    print(d_strut[1086])
    diff = abs(d_wing - d_strut)
    idx = np.argmin(diff)
    
    return F_strut_array[idx], diff[idx]


F_strut = np.arange(0, 1500000, 1000)
width = 1
optimum = indet_sys(F_strut, width, gamma, L_strut, A_S, A_E)
print(optimum)

#F_strut = np.arange(optimum - 2000, optimum + 2000, 0.1)
#optimum = deter_strut(F_strut, width, gamma, L_strut, A_S, A_E)
##print(optimum)
#
#F_str = optimum
#X = np.arange(0, 30, 1)
#
#d_lift = deter_d_distr(A_L, X, l_distr)
#d_weight = deter_d_distr(A_W, X, -w_distr)
#d_strut = deter_d_force(A_S, X, np.sin(gamma) * F_str)
#d_engine = deter_d_force(A_E, X, W_eng)
#d = d_lift + d_weight + d_strut + d_engine
#
#x_plot = np.arange(30, 0, -1)
#plt.plot(x_plot, d)
#plt.show()
#
#print(d[np.argmax(d)])