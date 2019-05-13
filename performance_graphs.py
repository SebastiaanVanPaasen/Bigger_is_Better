# -*- coding: utf-8 -*-
"""
Created on Mon May 13 16:00:12 2019

@author: nikki
"""
#-------------------------------MODULES---------------------------------------
import numpy as np
import matplotlib.pyplot as plt

#-------------------------------INPUTS-----------------------------------------
# S         =       Surface area 
#CD0        =       Zero lift drag coefficient
#A          =       Aspect ratio
#e          =       Oswald efficiency factor 
#R_des      =       Design range [m]
#Tcr        =       Thrust setting during cruise
#T_W        =       Thrust to weight ratio during cruise
#R_design   =       design range in m
#Wcr        =       Cruise weight in N
#F          =       Fuel flow during cruise kg/s


#------------------------------VERIFICATION DATA--------------------------------

"""Inputs unit test based on B737 MAX-8"""

MLW  = 69308.*9.81

Tcr = 117.3e03  #N
A = 8.45
e = 0.85
CD0 = 0.020
g = 9.81
S = 124.5 

F = 2000./3600. #Fuel flow kg/s

#assume weight in cruise is at least MLW to allow for a save landing
Wcr = MLW

#------------------------------DEFINITIONS-----------------------------------

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
        h0 = 11000.     #m
        rho = rho0*np.e**((-g*M*(h-h0))/(R*T))
        
    return rho
    

def drag_plot(h,S,A,e,Wcr):         #Drag VS velocity graph
    V = np.linspace(200,1200,800)   #V range in km/h
    CL_list = []
    CD_list = []
    for i in V:
        CL = Wcr / (0.5*ISA_density(h)*(i/3.6)**2*S)
        CD = CD0 + (CL**2 / (np.pi*A*e))
        CL_list.append(CL)
        CD_list.append(CD)
    
    D = []                  #drag at different airspeeds in kN
    for j in range(len(CL_list)):
        D.append((Wcr * (CD_list[j]/CL_list[j]))/1000.)
        
    #plt.plot(V, D)                   #V-D plot for a given altitude and varying speed
    #plt.xlabel("Speed [km/h]")
    #plt.show()
    
    k = D.index(min(D))            
    return V[k] ,D                    #speed at which the minimum drag is experienced
    

def power_plot(h,Wcr,S,Tcr):        #input weight in Newtons!!!
    V = np.linspace(200,1200,800)   #V range in km/h
    Pr_list = []                    #power required list in kW
    Pa_list = []                    #Power available list in kW
    for i in V:
        CL = Wcr / (0.5*ISA_density(h)*(i/3.6)**2*S)
        CD = CD0 + (CL**2 / (np.pi*A*e))
        Pr = Wcr*np.sqrt((2*Wcr*CD**2)/(S*ISA_density(h)*CL**3))
        
        Pr_list.append(Pr/1000.)
        Pa_list.append((i/3.6)*Tcr/1000.)
        
    
    Pc_list = []                    # Excess power used to climb in kW
    for i in range(len(Pa_list)):
        Pc = Pa_list[i] - Pr_list[i]
        Pc_list.append(Pc)
    
    RC_list = []                    # rate of climb is proportional to excess power
    for j in Pc_list:
        RC_list.append((j*1000)/Wcr)
    
    #Pr, Pa and RC plots for a given altitude
    #plt.plot(V,Pr_list,label = "P_req")
    #plt.plot(V,Pa_list, label = "P_av")
    #plt.legend(loc = "upper right")
    #plt.xlabel("Airspeed [km/h]")
    #plt.ylabel("Power [kW]")
    #plt.show()

    #plt.plot(V,RC_list)
    #plt.xlabel("Airspeed [km/h]")   
    #plt.ylabel("Rate of climb [m/s]")
    #plt.show()

    
    k = RC_list.index(max(RC_list))
    return V[k], max(RC_list),Pr_list,Pa_list,V  #peed at which the maxmimum RC is obtained 


#----------------------------MAIN PROGRAM-----------------------------------
"""Required power wrt altitud and speed at min Drag"""
dh = 100                #step size in altitude
H = range(7000,12000,dh)#altitude range

V_minD = []             #V at minimum drag
for h in H:
    V_minDi = drag_plot(h,S,A,e,Wcr)[0]
    V_minD.append(V_minDi)

plt.subplot(211)   
plt.plot(H,V_minD)
plt.xlabel("Altitude [m]")
plt.ylabel("Airspeed at D_min [km/h]")


for h in H:
    Pr = power_plot(h,Wcr,S,Tcr)[2]
    V = power_plot(h,Wcr,S,Tcr)[4]
    
    
    plt.subplot(212)
    plt.plot(V,Pr)
    plt.xlabel("Airspeed [km/h]")
    plt.ylabel("Power required [kW] at different H")

plt.show()        
    
"""RC max and corresponding airspeed and min time to climb to H wrt altitude"""
V_RCmax = []            #list with V at which RC is max at different altitudes
RC_max = []             #list with RC max values at different altitudes
tmin = [0]              #list with cummelative time to climb at RC max
W_fuel = []             #fuel consumed in kg

Wfi = 0                 #fuel consumption
dh = 100                #step size in altitude
H = range(7000,12000,dh)#altitude range

#RC max with respect to altitude
for h in H:
    Wi = Wcr - Wfi          #take into account the weight variation due to fuel
    RC_maxi = power_plot(h,Wi,S,Tcr)[1]
    
    RC_max.append(RC_maxi)

    Wfi +=  F*g*(dh/RC_maxi)
    W_fuel.append(Wfi/g)
    
#V at RC_max with respect to altitude
Wfi = 0                     #fuel consumption
for h in H:
    Wi = Wcr - Wfi          #take into account the weight variation due to fuel
    V_RCi = power_plot(h,Wi,S,Tcr)[0]

    V_RCmax.append(V_RCi)
    
    RC_maxi = power_plot(h,Wi,S,Tcr)[1]
    Wfi +=  F*g*(dh/RC_maxi)


#Minimum time to climb with respect to altitude
Wfi = 0                     #fuel consumption
for h in H:
    Wi = Wcr - Wfi          #take into account the weight variation due to fuel
    RCi = power_plot(h,Wi,S,Tcr)[1]
    RCii = power_plot(h+dh,Wi,S,Tcr)[1]
    RC_ave = (RCi+RCii)/2.
    
    tmin.append(tmin[-1]+(dh/RC_ave))

    Wfi +=  F*g*(dh/RC_ave)


tmin.remove(0)
plt.subplot(221)
plt.plot(H,tmin)
plt.ylabel("T min to climb [s]")
plt.xlabel("Altitude [m]")

plt.subplot(222)
plt.plot(H,RC_max)
plt.ylabel("RC_max [m/s]")
plt.xlabel("Altitude [m]")

plt.subplot(223)
plt.plot(H,V_RCmax)
plt.ylabel("V at RC_max [km/s]")
plt.xlabel("Altitude [m]")

plt.subplot(224)
plt.plot(H,W_fuel)
plt.xlabel("Altitude [m]")
plt.ylabel("Fuel consumption [kg]")


plt.show()






    
    
    

        
        
    



