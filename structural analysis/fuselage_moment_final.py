import numpy as np
from matplotlib import pyplot as plt

## Iput variables
X_fuse = 49.3
Weight_total = -264441.6*3.75
Wing_y = 1464076*2 
F_Str_y = -1198802*2
H_Tail_y = -910*2 
V_Tail_y = -16986.8322
Wing_pos = 0.3*5.53+19.11 
strut_pos = 0.3*5.53+19.11 
H_Tail_pos = 49.3#50.87
V_Tail_pos = 45.43

## list were everything is defined
dx = 1
        
X_fuse = np.arange(0,X_fuse+dx,dx)    
#X_tip = np.arange(((X_fuse[1]-X_fuse[0])/2.),(X_fuse[-1]+dx),dx)
#

### Shear forces
Weight = (Weight_total/len(X_tip)) * np.ones(len(X_tip))
#
##Lift_total = sum(Weight) + F_Str_y - H_Tail_y + V_Tail_y
##Lift = (Lift_total/len(X_tip)) * np.ones(len(X_tip))
#
#
##Vy_root = sum(Lift) - sum(Weight) - F_Str_y + H_Tail_y + V_Tail_y
#
#
#
### Moments
#Lift_mom = Lift * X_fuse
#Weight_mom = Weight * X_fuse
#Wing_mom = Wing_y * (X_fuse - Wing_pos)
#Strut_mom = F_Str_y * (X_fuse - strut_pos)
#H_Tail_mom = H_Tail_y * (X_fuse- H_Tail_pos)
#V_Tail_mom = V_Tail_y * (X_fuse - V_Tail_pos)     
#
#Mz_root = sum(Lift_mom) - sum(Weight_mom) + Wing_mom - Strut_mom + H_Tail_mom - V_Tail_mom
#
#Mz_dist = [None] * len(X_fuse)#np.zeros(len(X_fuse))
#Vy_dist = [None] * len(X_fuse)#np.zeros(len(X_fuse))
#
#
#Mz_dist[0] = Mz_root
#Vy_dist[0] = Vy_root
#
### Calculating shear forces per section
#for i in range(len(X_fuse)-1):
#    Vy_section = Vy_dist[i] - Lift[i] + Weight[i]
#    
#    if X_fuse[i] > (Wing_pos) and X_fuse[i - 1] < (Wing_pos):
#        Vy_section -= Wing_y
#        
#    if X_fuse[i] > (strut_pos) and X_fuse[i - 1] < (strut_pos):
#        Vy_section += F_Str_y
#        
##    if X_fuse[i] > (H_Tail_pos) and X_fuse[i - 1] < (H_Tail_pos):
##        Vy_section -= H_Tail_y
#        
#    if X_fuse[i] > (V_Tail_pos) and X_fuse[i - 1] < (V_Tail_pos):
#        Vy_section += V_Tail_y
#    
#    Vy_dist[i + 1] = Vy_section

## Calculating moments per section
#for i in range(len(X_fuse)-1):
#    Mz_dist[i + 1] = Mz_root - Vy_root * X_fuse[i+1] + sum(((Lift - Weight) * (X_fuse[i+1] - X_tip))[:i])
#    
#    if X_fuse[i] > (Wing_pos): #and X_root[i - 1] < (L_wing - A_E):
#        Mz_dist[i + 1] += Wing_y * (X_fuse[i+1]- (X_fuse - Wing_pos))
#        
#    if X_fuse[i+1] > (strut_pos): #and X_root[i - 1] < (L_wing - A_S_L[idx]):
#        Mz_dist[i + 1] -= F_Str_y *(X_fuse[i+1]- (X_fuse - strut_pos))
#        
#    if X_fuse[i+1] > (H_Tail_pos): #and X_root[i - 1] < (L_wing - A_E):
#        Mz_dist[i + 1] += H_Tail_y * (X_fuse[i+1]- (X_fuse - H_Tail_pos))
#
#    if X_fuse[i+1] > (V_Tail_pos): #and X_root[i - 1] < (L_wing - A_E):
#        Mz_dist[i + 1] -= V_Tail_y * (X_fuse[i+1]- (X_fuse - V_Tail_pos))


plt.plot(X_fuse, Vy_dist)
plt.show()
#print(max(Mz_dist))