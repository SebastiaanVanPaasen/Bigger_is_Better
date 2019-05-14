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
plt.subplot(212)
plt.plot(SER(CF,CT,F)[1],SER(CF,CT,F)[0])
plt.title("Economic specific range SER""")
plt.xlabel("Cruise speed [km/h]")
plt.ylabel("SER [km/N]""")

plt.show()

























