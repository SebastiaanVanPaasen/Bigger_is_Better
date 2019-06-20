# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 14:19:06 2019

@author: nvanluijk
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 15:36:13 2019

@author: nvanluijk
"""
import numpy as np
import matplotlib.pyplot as plt

#------------------------------INPUTS--------------------------------------------
m = 1.2

"""Design inputs"""
Wto = 1555182.652 
Wland = 1555182.652 - 0.8*171848.0064
Wcr = 1555182.652 - 0.4*171848.0064
T0 = 406.*1000

Vcr = 235.4151
S = 240.97
A = 15
e = 0.6
CD0 = 0.31
CD = 0.025

CL_maxcr = 1.62
CLmax_land = 3.21
CLmax_TO = 2.94
Vs_land = np.sqrt((2*Wland)/(1.225*S*CLmax_land))
Vs_TO = np.sqrt((2*Wto)/(1.225*S*CLmax_TO))

Va = 1.3*Vs_land
Vtd = 1.15*Vs_land
V_LOF = 1.2*Vs_TO 

"""B777-200 inputs"""
#Wcr = (242670. - 0.*52160)*9.81
#T0 = 342.5*1000*2
#V = 244.87
#S = 427.8
#A = 8.67
#e = 0.6
#CD0 = 0.021
#CD = 0.025

#----------------------------DEFINITIONS--------------------------------------      
""" ISA definitions""" 
def ISA_density(h):      # enter height in m
    M = 0.0289644       #kg/mol molar mass of Earth's air
    R = 8.3144590       #universal gas constant Nm/MolK
    g = 9.81
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
    
def ISA_temp(h):
    if h < 11000:
        T = 288.15 - 0.0065*h   #in Kelvin
        return T
    if h >= 11000:
        return 216.65    #in Kelvin
    
def Mach(V,h):                  #enter V in m/s
    gamma = 1.4                 #enter h in m
    R = 287 #J/kg/K
    a = np.sqrt(gamma*R*ISA_temp(h))
    M = V/a
    return M 

def Vel(M,h):
    gamma = 1.4                 #enter h in m
    R = 287 #J/kg/K
    a = np.sqrt(gamma*R*ISA_temp(h))
    V = a*M
    return V

#Varying thrust with velocity
def T_alt(T0,h):
    T = T0*(ISA_density(h)/ISA_density(0))**m
    return T


"""Climb performance"""
#Normal rate of climb at a given altitude and airspeed (steady, no acceleration)
def RC(W,T0,V,S,A,e,CD0,h):
    k = 1. / (np.pi*A*e)
    T = T_alt(T0,h)
    RC = V*((T/W) - 0.5*ISA_density(h)*V**2 * (S/W)*CD0 - (W/S)*((2*k)/(ISA_density(h)*V**2)))
    return RC

def RC_unsteady(W,T0,V,S,h,CD):  #at const. EAS thus accelerating during flight
    T = T_alt(T0,h)
    M = Mach(V,h)
    vg_dvdh = 0.5668*M**2            #Constant EAS in tropospere
    D = 0.5*ISA_density(h)*V**2*S*CD
    C = ((T-D)*V)/(W*(1. + vg_dvdh))
    return C


#Maximum climb angle (or steepest climb) and corresponding airspeed and RC 
def steep_climb(T0,W,S,CD0,A,e,h,V):
    T = T_alt(T0,h)
    k = 1. / (np.pi*A*e)
    theta_max = np.arcsin((T/W) - np.sqrt(4*CD0*(k))) 

    V_theta_max = np.sqrt((2/ISA_density(h))*(k/CD0)**0.5 * (W/S)*np.cos(theta_max))
    RC_max_theta = V_theta_max*np.sin(theta_max)
    
    return theta_max*(180/np.pi), V_theta_max,RC_max_theta   #return in degrees

#--------------------------------------MAIN PROGRAM------------------------------
"""Steepest climb"""
V = np.arange(50,300,5)              #Criose velocity in m/s
H = np.arange(5000,13000,1000) 

Vs_list = []
theta_max_list = []
V_thetamax_list = []
RC_maxtheta_list = []
Vmin_list = []

for h in H: 
    Vs = np.sqrt((2*Wcr)/(ISA_density(h)*S*CL_maxcr))             #to indentify whether V_theta_max is obtainable
    Vs_list.append(Mach(Vs,h))                                  # Since often theta_max cannot be reached since req. V is lower than Vstall

    theta_max = steep_climb(T0,Wcr,S,CD0,A,e,h,Vcr)[0]
    V_theta_max = steep_climb(T0,Wcr,S,CD0,A,e,h,Vcr)[1]
    RC_theta_max =  steep_climb(T0,Wcr,S,CD0,A,e,h,Vcr)[2] 

    theta_max_list.append(theta_max)
    V_thetamax_list.append(Mach(V_theta_max,h))
    RC_maxtheta_list.append(RC_theta_max)

plt.figure(5)   

plt.subplot(221)
plt.plot(H,theta_max_list)
plt.title("Max. climb angle")
plt.xlabel("Altitude [m]")
plt.ylabel("Climb angle [deg]")
plt.grid(True)

plt.subplot(222)
plt.plot(H,V_thetamax_list, label = "V at theta max")
plt.plot(H,Vs_list,"--",label = "Stall speed")
plt.title("Required speed at max. climb angle")
plt.xlabel("Altitude [m]")
plt.ylabel("Mach number")
plt.grid(True)
plt.legend()

plt.subplot(223)
plt.plot(H,RC_maxtheta_list)
plt.title("Rate of climb at max. climb angle")
plt.xlabel("Altitude [m]")
plt.ylabel("Rate of climb [m/s]")
plt.grid(True)

"""Steady Climb rate"""
#RC_app = RC(Wland,T0,Va,S,A,e,CD0,0)
#RC_land = RC(Wland,T0,Vtd,S,A,e,CD0,0)
#RC_to = RC(Wto,T0,V_LOF,S,A,e,CD0,0) #at const. EAS thus accelerating during flight
#
#theta_app = np.arcsin(RC_app/Va)*(180/np.pi)
#theta_land = np.arcsin(RC_land/Vtd)*(180/np.pi)
#theta_to = np.arcsin(RC_to/V_LOF)*(180/np.pi)

#theta_app = RC_app/Va*100
#theta_land = RC_land/Vtd*100
#theta_to = RC_to/V_LOF*100

#print (theta_app,theta_land,theta_to)
#H = np.arange(1000,11000,1000)
#V = np.arange(80,350,20)
#
#plt.figure(1)
#RC_max = [] 
#M_RC_max = []
#V_RC_max = []
#H_RC_max = []
#
#for h in H:
#    RC_list = []
#    M_list = []
#    V_list = []
#    
#    for v in V:
#        RC_list.append(RC(Wcr,T0,v,S,A,e,CD0,h))
#        M_list.append(Mach(v,h))
#        V_list.append(v)
#        
#    plt.plot(M_list,RC_list, label = "%s m" %h)
#    
#    k = RC_list.index(max(RC_list))
#    RC_max.append(max(RC_list))
#    V_RC_max.append(V_list[k])
#    M_RC_max.append(M_list[k])
#    H_RC_max.append(h)
#    
#
#plt.plot(M_RC_max,RC_max," ko",label = " Max. RC" )
##plt.title("Steady rate of climb" )
#plt.xlabel("Mach number", fontsize = "x-large" )
#plt.ylabel("Rate of climb [m/s]", fontsize = "x-large" )
#plt.legend(fontsize = "x-large" )  
#plt.grid(True)
##plt.show()
#
#plt.figure(2)
#plt.plot(H,RC_max)
##plt.hlines(10.16,H[0],H[-1])
##plt.title("Max. steady rate of climb" )
#plt.xlabel("Altitude [m]",fontsize = "x-large"  )
#plt.ylabel("Maximum rate of climb [m/s]" , fontsize = "x-large" )
#plt.grid(True)

#plt.figure(3)
#plt.plot(H,M_RC_max)
##plt.title("Mach number at Max. steady RC" )
#plt.xlabel("Altitude [m]")
#plt.ylabel("Mach number")
#plt.grid(True)
    
#plt.show()    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




