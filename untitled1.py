import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from class_I.lift_distr import get_correct_data, lift_distribution, make_avl_file

D_fus = 7.3

A_strut = 0.25 * np.pi * ((5 / 1000) ** 2)
E_strut = 181 * (10 ** 9)  

E_wing = 69 * (10 ** 9)
I_wing = 0.000499
L_wing = 30

cr = 6.141
ct = 1.826

rho_cr = 
V_cr = 

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


def deter_d_distr(applic, x, distr):
    Y_distr = (distr / (24 * E_wing * I_wing)) * ((L_wing - applic) ** 3) * (3 * L_wing + applic)
    Theta_distr = (-distr / (6 * E_wing * I_wing)) * ((L_wing - applic) ** 3) * x
    Mac = np.array([])
    for i in range(len(Theta_distr)):
        if X[i] < applic:
            Mac = np.append(Mac, 0)
        else:
            mac = (distr / (24 * E_wing * I_wing)) * ((x[i] - applic) ** 4)
            Mac = np.append(Mac, mac)
        
    d_distr= Y_distr+ Theta_distr + Mac
    
    return d_distr


def deter_d_force(applic, x, force):
    Y_force = (-force/ (6 * E_wing * I_wing)) * (2 * (L_wing ** 3) - 3 * (L_wing ** 2) * applic + applic ** 3)
    Theta_force= (force / (2 * E_wing * I_wing)) * ((L_wing - applic) ** 2) * x
    Mac = np.array([])
    for i in range(len(Theta_force)):
        if X[i] < applic:
            Mac = np.append(Mac, 0)
        else:
            mac = (-force / (6 * E_wing * I_wing)) * ((x[i] - applic) ** 3)
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
        
    
#    def calc_volume(self):
#        self.volume = self.width * self.chord * (self.tc * self.chord)
        

def deter_lift(cl_curve, x, dx):
    CL = np.aray([])
    S = np.array([])
    
    for i in range(len(x)):
        CL = np.append(CL, cl_curve(x[i]))
        S = np.append(S, calc_area(x[i], dx))
        
        
    qs = 0.5 * rho_cr * (V_cr ** 2) * S
    L = qs * CL
    print(L)
    
    return L


cl = 0.37
l_distr = 25295.5
w_distr = 3677
W_eng = 73144
force = l_distr, w_distr, W_eng



A_L = 0
A_W = 0
A_E = 23
A_S = 5
applic = A_L, A_W, A_E, A_S

cl_curve, cd_curve = get_data(cl)


L_strut = np.sqrt(D_fus ** 2 + (L_wing - A_S) ** 2)
gamma = np.arctan(D_fus / (L_wing - A_S))
#print(gamma)

F_strut = np.arange(0, 1500000, 1000)
optimum = deter_strut(F_strut, gamma, L_strut, applic, force)
#print(optimum)

F_strut = np.arange(optimum - 2000, optimum + 2000, 0.1)
optimum = deter_strut(F_strut, gamma, L_strut, applic, force)
#print(optimum)

F_str = optimum
X = np.arange(0, 30, 1)

d_lift = deter_d_distr(A_L, X, l_distr)
d_weight = deter_d_distr(A_W, X, -w_distr)
d_strut = deter_d_force(A_S, X, np.sin(gamma) * F_str)
d_engine = deter_d_force(A_E, X, W_eng)
d = d_lift + d_weight + d_strut + d_engine

x_plot = np.arange(30, 0, -1)
plt.plot(x_plot, d)
plt.show()

print(d[np.argmax(d)])