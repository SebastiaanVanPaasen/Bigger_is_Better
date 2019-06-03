# -*- coding: utf-8 -*-
"""
Created on Tue May 14 11:13:37 2019

@author: nikki
"""

#-------------------------------MODULES---------------------------------------
import numpy as np
import matplotlib.pyplot as plt

#-------------------------------INPUTS-----------------------------------------
"""Aicraft configuration"""
#R_des      =       Design range [m]
#Tcr        =       Thrust during cruise [N]

"""Economic factors"""
#CT         =       Cost of time in $/h
#CF         =       Cost of fuel in $/N


#------------------------STATISTICAL INPUTS----------------------------------

Ct           = 12e-06      #Specific fuel conspumtion [kg/N/s] from B737MAX 


#------------------------------VERIFICATION DATA--------------------------------

"""Inputs unit test based on B737 MAX-8"""
R_des = 7000 #[km]
Tcr = 117.3e03  #N
CT = 200 #$/hr
CF  = 72.23 #$/N
F = Ct*Tcr*3600.*9.81 #N/hr
#-----------------------------------DEFINITIONS-------------------------------
"""Block speed and transport productivity"""

def ISA_density(h):      # enter height in m
    M = 0.0289644       #kg/mol molar mass of Earth's air
    R = 8.3144590       #universal gas constant Nm/MolK
    
    if h < 11000:
        rho0 = 1.225   #kg/m^3
        T = 288.15     #K
        h0 = 0.         #m
        a = -0.0065    #K/m
        rho = rho0 * (T/(T + a*(h-h0)))**(1. + ((g*M)/(R*a)))
        
    if h >= 11000:
        rho0 = 0.36391 #kg/m^3
        T = 216.65     #K
        h0 = 11000.    #m
        rho = rho0*np.e**((-g*M*(h-h0))/(R*T))
        
    return rho
    
def ISA_temp(h):
    if h < 11000:
        T = 288.15 - 0.0065*h   #in Kelvin
        return T
    if h >= 11000:
        return 216.65    #in Kelvin
        
def Mach(V,h):           #enter V in km/h and h in meters
    gamma = 1.4
    R = 287 #J/kg/K
    a = np.sqrt(gamma*R*ISA_temp(h))
    M = (V/3.6)/a
    return M 


def block_speed(Vcr):                 #enter in km/s
    R = range(500,10000,100)         #stage length
    dt = 50*60.                     #min (time between starting and shutting down engines, no ground operations)
    
    V_block = []                    #result in m/s
    for r in R:
        V_B = (r*1000) / (((r*1000)/(Vcr/3.6)) + dt)
        V_block.append(V_B)
    
    return V_block, R


"""SER: economic specific range"""
def SER(CF,CT,F):
    CI = CT/CF
    V = range(200,1000,100)
    SER = []
    for v in V:
        SER.append((v)/(CI+F))

    return SER, V
    
    

#-------------------------------MAIN PROGRAM---------------------------------
V = range(200,1000,100)
R = range(500,2000,100) 

#Block speed diagram
for v in V:
    plt.subplot(211)
    plt.plot(block_speed(v)[1],block_speed(v)[0],label='%s cruise speed [km/h]' % v)
    plt.title('Block speed')
    plt.xlabel("Range [km]")
    plt.ylabel("Speed [km]/h")
    
plt.legend(loc = 'upper left') 

#SER diagram
dh = 500                #step size in altitude
H = range(7000,12000,dh)#altitude range

for h in H:   
    SER_lst = SER(CF,CT,F)[0]
    M = []
    for i in range(len(SER(CF,CT,F)[1])):
        Mi = Mach(SER(CF,CT,F)[1][i],h)
        M.append(Mi)
    
    plt.subplot(212)
    plt.plot(M,SER_lst,label='%s altitude [m]' % h)
    plt.title("Economic specific range SER""")
    plt.xlabel("Cruise Mach number")
    plt.ylabel("SER [km/N]""")

plt.show()

























