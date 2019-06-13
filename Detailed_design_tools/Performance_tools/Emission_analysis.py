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
#Tc = total inlet combustor temp. [K]
#R_des = design range [km]

SAR = 0.0094
n_eng = 2.
n_pax = 200.
hcr = 10000.
Mcr = 0.79
MTOW = 78245.
R_des = 4500. 

Tc = 711.5    #from cycle calculation

#------------------------------------DEFINITIONS-------------------------------
def ISA_temp(h):
    if h < 11000:
        T = 288.15 - 0.0065*h   #in Kelvin
        return T
    if h >= 11000:
        return 216.65    #in Kelvin
        
def ISA_press(h):
    R = 8.3144590               #universal gas constant Nm/MolK
    M = 0.0289644               #kg/mol molar mass of Earth's air
    g = 9.80665                 #gravitationa constant 
    
    p0 = 101325.0
    T0 = 288.15 #[K]
    a = -0.0065 
    T1 = ISA_temp(11000)
    
    if h == 0:
        p = p0
    elif h <= 11000:
        p = p0*(T0/(T0 + a*h))**((g*M)/(R*a))
        #p = p0*(T/T0)**(-g/(a*R))
    else:
        p = 22632.10 *np.exp((-g*M*(h - 11000.))/(R*T1))
        
    return p
    

def Vel(M,h):
    gamma = 1.4                 #enter h in m
    R = 287 #J/kg/K
    a = np.sqrt(gamma*R*ISA_temp(h))
    V = a*M
    return V

def LTO_emis(n_eng,Ffuel,TIM,EF):    #Ffuel is fuel flow in kg/s
    emission = n_eng*Ffuel*TIM*EF    #EF is emission factor in kg/kg
    return emission

def EF_NOx(Tc,h): #in kg/kg
    #NOx varies with altitude, temp. in cc
    EF_NOx = (10**(1. + 0.0032*(Tc - 581.25)))*np.sqrt(ISA_press(h)/ISA_press(0))
    return EF_NOx/1000.
    
def Cruise_emis(SAR,EF):#returns emission in kg/km-pax
    emission = SAR*EF
    return emission


#-----------------------------------LTO CYCLE-----------------------------------
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
NOx = EF_NOx(Tc,0.)
NMVOC = 0.00743 
CO = 0.02367
PM10 = 0.000064 

EF = [H2O,CO2,SO2,CH4,N2O,NOx,NMVOC,CO,PM10]
EF_comp = ['H2O','CO2','SO2','CH4','N2O','NOx','NMVOC','CO','PM10']

"""Tyre and brake wear emission factors: TNO source"""
PM10_tyre = 2.23e-07   #kg/kg MTOW
PM10_brake = 2.53e-07  #kg/kg MTOW


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
print()





#--------------------------------CRUISE----------------------------------------
"""Cruise emissions"""
#Known that CH4 emission in cruise is practically zero, we can neglect it
#Consider the main emissions during cruise: NOx, CO, VOC, SO2, CO2, H2O
#Note CO2 and H2O values are simliar for the ground operations
#While NOx varies with altitude and combustor temperature
#As values in cruise EF are highly uncertain, only the main emissions were 
#considere for the cruise phase
#As cruising time differs each missions, the emissions are calculated using the SAR
#This resulted in the emissions in kg/(km-pax)
"""Global warming potentials"""
#Source 2
GWP_CO2 = 1.
GWP_CH4 = 25.
GWP_N2O = 298.


"""Emission factors"""
#Source 1 [kg/kg fuel]
H2O = 1.237 
#Source 3  [kg/kg fuel]
CO2 =3.15

#NOx, see ruijgrok elements of a/c pollution book
NOx = EF_NOx(Tc,hcr)


EF_comp = ['H2O','CO2','NOx','CO','NMVOC','NOx','SO2']

if hcr <= 9100.:
    EF = [H2O,CO2,NOx,0.0052,0.0009,0.0107,0.001]
elif hcr > 9100.:
    EF = [H2O,CO2,NOx, 0.005,0.0008,0.0108,0.001]


"""Emission calculations"""
print ("Cruise emissions [kg/km-pax] & kg")
Emissions_cruise = []

for i in range(len(EF)):
    emissioni = Cruise_emis(SAR,EF[i])
    Emissions_cruise.append(emissioni)
    
    #In case the design range is the entire cruise distance    
    total_emis_cr = emissioni*R_des*n_pax 
        
    print (EF_comp[i],emissioni,total_emis_cr)
    
 
"""Sensitivity EF of NOx"""
Tc_list = np.arange(500,1200,50)
H = np.arange(0,12500,500)

NOx_consth = []
for t in Tc_list: #for constant cruise altitude
    NOx = EF_NOx(t,hcr)
    NOx_consth.append(NOx)
    
NOx_constT = []
for h in H:
    NOx = EF_NOx(Tc,h)
    NOx_constT.append(NOx)
    
plt.subplot(121)
plt.plot(Tc_list, NOx_consth)
plt.title("EF NOx for const. h, varying Tc")
plt.grid(True)
plt.xlabel("Tc [K]")
plt.ylabel("EF NOx [kg/kg]")

plt.subplot(122)
plt.plot(H, NOx_constT)
plt.title("EF NOx for const. Tc, varying h")
plt.grid(True)
plt.xlabel("h [m]")
plt.ylabel("EF NOx [kg/kg]")
plt.show()  


#---------------------------------------RADIATIVE FORCING-------------------------
#To take into account the effect of the emissions during cruise as they do not directly
#influence the quality of life on the ground, while the LTO emissions do
#To take into account the phenomena which occur high in the atmosphere due to
#aircraft emissions -> contrails and cirrus clouds   
    
    

    





















