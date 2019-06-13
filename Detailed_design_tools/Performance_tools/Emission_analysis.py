# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 17:11:54 2019

@author: nikki
"""

import numpy as np
import matplotlib.pyplot as plt

#----------------------------------INPUTS--------------------------------------
#SAR = cruise fuel consumption in kg/km-pax
#n_eng = number of engines
#n_pax = number of passengers
#hcr = cruise altitude [m]
#Mcr = cruise Mach number
#MTOW = maximum tak-off weight in [kg]

SAR = 0.0094
n_eng = 2.
n_pax = 200.
hcr = 10000.
Mcr = 0.79
MTOW = 78245.


#-----------------------------------STANDARD INPUTS---------------------------
"""Time in mode [s]"""
TIM_TO = 0.7*60. 
TIM_climb = 2.2*60.
TIM_app = 4.*60.
TIM_idle = 26.*60.  #Includes taxiing and holding position

TIM = [TIM_TO,TIM_climb,TIM_app,TIM_idle]

"""Thrust settings as fractions"""
TS_TO = 1.
TS_climb = 0.85
TS_app = 0.3
TS_idle = 0.07
TS_cruise = 0.6

TS = [TS_TO,TS_climb,TS_app,TS_idle]
phase = ['TO','climb',"app","idle"]


"""Emission factors"""
#Source 1 [kg/kg fuel]
H2O = 1.237 

#Source 3 in [kg/kg fuel]
CO2 =3.15
SO2 = 0.000274 
CH4 = 0.00035 
N2O = 0.0001  
NOx = 0.01351 
NMVOC = 0.00743 
CO = 0.02367
PM10 = 0.000064 

EF = [H2O,CO2,SO2,CH4,N2O,NOx,NMVOC,CO,PM10]
EF_comp = ['H2O','CO2','SO2','CH4','N2O','NOx','NMVOC','CO','PM10']

"""Tyre and brake wear emission factors: TNO source"""
PM10_tyre = 2.23e-07   #kg/kg MTOW
PM10_brake = 2.53e-07  #kg/kg MTOW


"""Global warming potentials"""
#Source 2
GWP_CO2 = 1.
GWP_CH4 = 25.
GWP_N2O = 298.





#------------------------------------DEFINITIONS-------------------------------
def ISA_temp(h):
    if h < 11000:
        T = 288.15 - 0.0065*h   #in Kelvin
        return T
    if h >= 11000:
        return 216.65    #in Kelvin

def Vel(M,h):
    gamma = 1.4                 #enter h in m
    R = 287 #J/kg/K
    a = np.sqrt(gamma*R*ISA_temp(h))
    V = a*M
    return V

def LTO_emis(n_eng,Ffuel,TIM,EF):    #Ffuel is fuel flow in kg/s
    emission = n_eng*Ffuel*TIM*EF    #EF is emission factor in kg/kg
    return emission


#-----------------------------------MAIN PROGRAM--------------------------------
"""Constant variables calculation"""
#Convert SAR to fuel flow in kg/s for cruise conditions so at TS_cruise
Vcr = Vel(Mcr,hcr)                 #m/s
Ffuel_cr = (SAR/1000.)*n_pax*Vcr   

#Convert Ffuel for a TS of 100
Ffuel_max = (Ffuel_cr/(TS_cruise*100.))*100.


"""LTO cycle emissions"""
Emis_tot = []

for i in range(len(phase)):
    E_list = []
    for j in range(len(EF)):
        EF_emis = EF[j]
        TIM_phase = TIM[i]
        TS_phase = TS[i]
        
        Ffuel = Ffuel_max*TS_phase
        emission = LTO_emis(n_eng,Ffuel,TIM_phase,EF_emis)        
        E_list.append(emission)

    Emis_tot.append(E_list)
    
    
print ("Order of emissions ([kg]):",EF_comp)
print ()
    
for k in range(len(Emis_tot)):
    print ("Phase= ",phase[k])
    print (Emis_tot[k])
    print ()
    

"""Tyre and break wear additional emissions"""
Tyre_wear = PM10_tyre*MTOW
Brake_wear = PM10_brake*MTOW

print ("PM10 emissions due to brake and tyre wear:", Tyre_wear, Brake_wear)


"""Cruise emissions"""





















