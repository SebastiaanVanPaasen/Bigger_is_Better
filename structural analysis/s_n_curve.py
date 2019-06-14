import numpy as np
from math import *
from matplotlib import pyplot as plt
##Fatigue 

##Aluminium properties:
S_ult =(0.667)*(570*10**6)
material = 'Aluminium'
load_type = 'bending'
Kt = 2
##Properties for low cycles

## define materials properties witn boundary conditions!
if material == 'Steel':
    if S_ult<(1400*10**6):
        Se = 0.5*S_ult
        Sf = 'NvT'
    elif S_ult>=(1400*10**6):
        Se = 700*(10**6)
        Sf = 'NvT'
elif material == 'Irons':
    if S_ult<(400*10**6):
        Se = 0.4*S_ult
        Sf = 'NvT'
    elif S_ult>=(400*10**6):
        Se = 160*(10**6)
        Sf = 'NvT'
elif material == 'Aluminium':
    if S_ult<(330*10**6):
        Sf = (0.4*S_ult)/Kt
        Se = 'NvT'
    elif S_ult>=(330*10**6):
        Sf = (130*10**6)/Kt
        Se = 'NvT'
elif material == 'Copper alloys':
    if S_ult<(280*10**6):
        Sf = (0.4*S_ult)/Kt
        se = 'NvT'
    elif S_ult>=(280*10**6):
        Sf = (100*10**6)/Kt
        Se = 'NvT'



if Sf == 'NvT':
    Nf = np.arange((10**3),(10**9),1000)
    Sn_complete = []
    N_Se = (1*10**6)
    if load_type == 'bending':
        Sm = 0.9*S_ult
        N_Sm = (1*10**3)
    elif load_type == 'axial load':
        Sm = 0.75*S_ult
        N_Sm = (1*10**3)
    z = -3
    b = (1/z)*log10(Sm/Se)
    a = 10**(log10(Sm)-b*log10(N_Sm))
    Sn_complete.append(Sm)
    for i in range(1,len(Nf)):
        if Nf[i]<=N_Se:
            Sn = a*(Nf[i]**b)
            Sn_complete.append(Sn)
        elif Nf[i] > N_Se:
            Sn_complete.append(Sn)

elif Se == 'NvT':
    Nf = np.arange((10**3),(10**9),1000)
    Sn_complete = []
    N_Sf = 5*10**8
    if load_type == 'bending':
        Sm = 0.9*S_ult
        N_Sm = (1*10**3)
    elif load_type == 'axial load':
        Sm = 0.75*S_ult
        N_Sm = (1*10**3)
    z = -5.699
    b = (1/z)*log10(Sm/Sf)
    a = 10**(log10(Sm)-b*log10(N_Sm))
    Sn_complete.append(Sm)
    for i in range(1,len(Nf)):
        Sn = a*(Nf[i]**b)
        Sn_complete.append(Sn)

life_time = (round((6*7*365*27*8)/1000)*1000)
Nf = list(Nf)
pos = Nf.index(life_time)
stress_life_time = Sn_complete[pos]

cycles = [life_time,life_time,0]
load = [0,stress_life_time,stress_life_time]


            
plt.plot(Nf,Sn_complete,label='Aluminium S-N Curve')
plt.plot(cycles,load, linestyle='dashed',label='Fatigue stress at design range')
plt.ylim(0,(400*10**6))
plt.legend(loc='upper right')
plt.xscale('log')
plt.xlabel('Number of cycles')
plt.ylabel('Stress')
plt.title('S-N Curve')
plt.show()      
        
        
        