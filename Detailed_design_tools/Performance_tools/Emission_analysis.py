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

"""Design"""
SAR = 0.0094
n_eng = 2.
n_pax = 200.
hcr = 11000.
Mcr = 0.79
MTOW = 82919.
R_des = 1100. 

Tc = 711.5    #from cycle calculation 

APU = True

"""Ref. aircraft"""
hcr_ref = 10000.

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
    
def Cruise_emis(SAR,EF):#returns emission in kg/km
    emission = SAR*EF*n_pax
    return emission
    
def RF_contrails(RF,h):
    RF_red = RF - 0.002935*((36000*0.3048)-h)
    return RF_red

def RF_NOx(RF,h):
    if h < 9144.:
        k = RF
    elif 9144. < h <= 9296.4:
        k = RF*1.19
    elif 9296.4 < h <= 9906.:
        k = RF*1.38
    elif 9906. < h <= 10515.6:
        k = RF*1.76
    elif 10515.6 < h:
        k = RF*2.
    return k    


#-----------------------------------LTO CYCLE-----------------------------------
"""Time in mode [s]"""
TIM_TO = 0.7*60. 
TIM_climb = 2.2*60.
TIM_app = 4.*60.
TIM_idle = 26.*60.  #Includes taxiing and holding position
TIM_APU = 0. 

if APU == True:
    TIM_APU = 30.*60.

TIM = [TIM_TO,TIM_climb,TIM_app,TIM_idle,TIM_APU]

"""Thrust settings as fractions"""
TS_TO = 1.
TS_climb = 0.85
TS_app = 0.3
TS_idle = 0.07
TS_cruise = 0.6

TS = [TS_TO,TS_climb,TS_app,TS_idle]
phase = ['TO','climb',"app","idle","APU"]


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

#Fuel flow of APU
Ffuel_APU = 2./60. #kg/s

"""LTO cycle emissions"""
Emis_tot = []

for i in range(len(phase)):
    E_list = []
    for j in range(len(EF)):
        if i <= 3.:
            EF_emis = EF[j]
            TIM_phase = TIM[i]
            TS_phase = TS[i]
            
            Ffuel = Ffuel_max*TS_phase
            emission = LTO_emis(n_eng,Ffuel,TIM_phase,EF_emis) 
            
            E_list.append(emission)
        elif i == 4.:
            EF_emis = EF[j]
            TIM_phase = TIM[i]
                
            emission = EF_emis*TIM_APU*Ffuel_APU
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

"""Emission factors"""
#Source 1 [kg/kg fuel]
H2O = 1.237 
#Source 3  [kg/kg fuel]
CO2 =3.15
N2O = 0.0001 

#NOx, see ruijgrok elements of a/c pollution book
NOx = EF_NOx(Tc,hcr)

EF_comp = ['H2O','CO2','NOx','CO','NMVOC','SO2',"N2O"]

EF = [H2O,CO2,NOx,0.005,0.0007,0.001,1e-04]





"""Emission calculations"""
print ("Cruise emissions [kg/km] & kg")
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
#Only consider the driver green house gases, like H2O, CO2 and NOx, which have the
#largest impact on the ozone layer

"""RF values in mW/m^2"""
#RF for each emission : best estimate [0], low [1], high [2], 90% confidence interval
CO2 = [28.,15.2,40.8]
CH4 = [0.,0.,0.]
NOx = [12.6,3.8,15.7]
H2O = [2.8,0.39,20.3]
contrails = [11.8,5.4,25.6]
AIC = [33.-contrails[0],12.5-contrails[1],86.7-contrails[2]]
SO4 = [-4.8,-0.79,29.3] #occur due to chemical reactions with SO2
soot = [3.4,0.56,20.7]
    
H = np.arange(5000,11500,500)

RF_tot_mean = []
RF_tot_low = []
RF_tot_high = []

for h in H:
    RF_tot_mean.append(CO2[0]+CH4[0]+RF_NOx(NOx[0],h)+H2O[0]+RF_contrails(contrails[0],h)+AIC[0]+SO4[0]+soot[0])
    RF_tot_low.append(CO2[1]+CH4[1]+RF_NOx(NOx[1],h)+H2O[1]+RF_contrails(contrails[1],h)+AIC[1]+SO4[1]+soot[1])    
    RF_tot_high.append(CO2[2]+CH4[2]+RF_NOx(NOx[2],h)+H2O[2]+RF_contrails(contrails[2],h)+AIC[2]+SO4[2]+soot[2])
    
#RF of aircraft
RF_ac = CO2[0]+CH4[0]+RF_NOx(NOx[0],hcr)+H2O[0]+RF_contrails(contrails[0],hcr)+AIC[0]+SO4[0]+soot[0]
RF_ref = CO2[0]+CH4[0]+RF_NOx(NOx[0],hcr_ref)+H2O[0]+RF_contrails(contrails[0],hcr_ref)+AIC[0]+SO4[0]+soot[0]

RF_reduction = (RF_ac - RF_ref)/RF_ref * 100.
print ("RF reduction wrt ref. aircraft",RF_reduction)



plt.plot(hcr_ref,RF_ref,"o",label = "Reference aicraft")
plt.plot(hcr,RF_ac,'o',label = "Aircraft")
plt.plot(H,RF_tot_mean,label = "Best estimate")
#plt.plot(H,RF_tot_low, label = "Lower bound values")
#plt.plot(H,RF_tot_high,label = "Upper bound values")
plt.title("RF dependency on altitude")
plt.xlabel("Altitude [m]")
plt.ylabel("RF [mW/m^2]")
plt.grid(True)

plt.show()


    
    





















