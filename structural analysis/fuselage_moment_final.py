import numpy as np
from matplotlib import pyplot as plt

## Iput variables
L_f = 49.3
L_cabin = 35.2
X_cabin = 4

W_TO = 3.75 * 1520276
W_PL = 3.75 * 472485
W_wing = 156375
W_fuel = 170698

W_nac = 25875
W_eng = 76148
W_propsys = W_nac + W_eng

W_distr = W_TO - W_PL - W_wing - W_fuel - W_propsys


# Strut characteristics

Vy_wing = 2 * 1464076
Vy_strut = 2 * 1181491

Wing_pos = 0.3*5.53+19.11 
strut_pos = 0.3*5.53+19.11 

# Tail characteristics

H_Tail_pos = 49.3  # 50.87
V_Tail_pos = 45.43

Vy_tail = Vy_wing + Vy_strut - W_PL - W_distr

## list were everything is defined
dx = 0.1 
X_fuse = np.arange(0, L_f + dx, dx)    
X_tip = np.arange(L_f - dx / 2, 0, -dx)

w_distr = np.ones(len(X_fuse) - 1) * W_distr / (L_f / dx)
#w_distr = W_distr / (L_f / dx)
w_pl = np.zeros(len(X_fuse) - 1)
w_pl_d = W_PL / (L_cabin / dx)


Vy_dist = np.zeros(len(X_fuse))
Mz_dist = np.zeros(len(X_fuse))

Vy_dist[0] = 0
Mz_dist[0] = 0

w_pl_tot = 0
w_distr_tot = 0

## Calculating shear forces per section
for i in range(len(X_fuse) - 1):
    Vy_section = Vy_dist[i] - w_distr[i]
    w_distr_tot += w_distr[i]
    
    if X_fuse[i] > Wing_pos and X_fuse[i - 1] < Wing_pos:
        Vy_section += Vy_wing
        
    if X_fuse[i] > strut_pos and X_fuse[i - 1] < strut_pos:
        Vy_section += Vy_strut
        
    if X_fuse[i] >= X_cabin and X_fuse[i + 1] <= (X_cabin + L_cabin):
        w_pl[i] = w_pl_d
        Vy_section -= w_pl[i]
        w_pl_tot += w_pl[i]
    else:
        w_pl[i] = 0
        
    if X_fuse[i + 1] >= H_Tail_pos and X_fuse[i] <= H_Tail_pos:
        Vy_section -= Vy_tail

    
    Vy_dist[i + 1] = Vy_section

## Calculating moments per section
for i in range(len(X_fuse)- 1):
    Mz_dist[i + 1] = sum(((w_distr + w_pl) * (X_fuse[i+1] - X_tip[::-1]))[:i])
    
    if X_fuse[i+1] > Wing_pos: #and X_root[i - 1] < (L_wing - A_E):
        Mz_dist[i + 1] -= Vy_wing * (X_fuse[i+1] - Wing_pos)
        
        
    if X_fuse[i+1] > strut_pos: #and X_root[i - 1] < (L_wing - A_S_L[idx]):
        Mz_dist[i + 1] -= Vy_strut *(X_fuse[i+1] - strut_pos)
        
    if X_fuse[i+1] > H_Tail_pos: #and X_root[i - 1] < (L_wing - A_E):
        Mz_dist[i + 1] += Vy_tail * (X_fuse[i+1] - H_Tail_pos)
        
Mx_wing = Mz_dist[-1] * (Vy_wing / (Vy_tail + Vy_wing))
Mx_tail = Mz_dist[-1] * (Vy_tail / (Vy_tail + Vy_wing))

for i in range(len(X_fuse) - 1):
    
    if X_fuse[i + 1] > Wing_pos:
        Mz_dist[i + 1] -= Mx_wing
        
    if X_fuse[i + 1] > H_Tail_pos:
        Mz_dist[i + 1] -= Mx_tail



plt.plot(X_fuse, Mz_dist)
plt.show()
#print(max(Mz_dist))