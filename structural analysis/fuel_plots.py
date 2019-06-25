import numpy as np
import constants_and_conversions as cc
import matplotlib.pyplot as plt


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

def calc_chord(x):
    return cr - ((cr - ct) / L_wing) * x

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

dx = 0.1
X_root = np.arange(0, L_wing+dx, dx) 
X_tip = np.arange(L_wing - dx / 2, 0, -dx)
x_start = [0, L_wing * 0.28, L_wing * 0.35]
weights = []

for i in range(len(x_start)):
    weights_v, Volumes = deter_weight(W_wing, X_tip[::-1], dx)  #
    fuel_weights_v = deter_fuel(W_fuel / 2, Volumes, Rho_fuel, X_tip[::-1], x_start[i])
    weights.append(fuel_weights_v)
    
plt.rcParams.update({"font.size": 20})

plt.subplot(1, 3, 1)
plt.plot(X_tip[::-1], weights[0])
plt.xlabel("X-position [m]")
plt.ylabel("Fuel weight [N]")
plt.title("Fuel distribution from root")

plt.subplot(1, 3, 2)
plt.plot(X_tip[::-1], weights[1])
plt.xlabel("X-position [m]")
plt.ylabel("Fuel weight [N]")
plt.title("Fuel distribution from 0.28 b/2")

plt.subplot(1, 3, 3)
plt.plot(X_tip[::-1], weights[2])
plt.xlabel("X-position [m]")
plt.ylabel("Fuel weight [N]")
plt.title("Fuel distribution from 0.35 b/2")

plt.show()