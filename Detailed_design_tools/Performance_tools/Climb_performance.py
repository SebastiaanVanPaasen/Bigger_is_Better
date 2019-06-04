# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 15:36:13 2019

@author: nvanluijk
"""
import numpy as np
import matplotlib.pyplot as plt

#------------------------------INPUTS--------------------------------------------
MTOW       =       78220*9.81
MFW        =       31594*9.81*0.75
Wcr        =       MTOW-0.3*MFW 
T0         =       127.62*1000*2  #@ SEA LEVEL!!!

A          =       9.44       
e          =       0.85 
S          =       124.60 
CD0        =       0.020
CL_max     =       2.2


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
    M = (V)/a
    return M 

""" Max and min velocities""" 
""" Maybe add thrust variation definition here too""" 
def Vmax(Wcr,Tmax,S,A,e,CD0,h):
    #In order to minimise mistakes, split up eq. in more variables
    k = 1. / (np.pi*A*e)
    Vmax = (((Tmax/Wcr)*(Wcr/S) + (Wcr/S)*np.sqrt((Tmax/Wcr)**2 - (4*CD0*k)))/(ISA_density(h)*CD0))**0.5
    return Vmax

def Vmin(Wcr,Tmax,S,A,e,CD0,h,CL_max):
    #In order to minimise mistakes, split up eq. in more variables
    k = 1. / (np.pi*A*e)
    Vmin = (((Tmax/Wcr)*(Wcr/S) - (Wcr/S)*np.sqrt((Tmax/Wcr)**2 - (4*CD0*k)))/(ISA_density(h)*CD0))**0.5
    
    Vs = np.sqrt((2*Wcr)/(ISA_density(h)*S*CL_max))
    
    if Vs > Vmin:
        return Vs
    if Vmin > Vs:
        return Vmin


""" Thrust required""" 
def Treq(W,V,S,A,e,h,CD0):
    CL = W / (0.5*ISA_density(h)*(V)**2*S)
    CD = CD0 + (CL**2 / (np.pi*A*e))
    
    Pr = np.sqrt((2*W**3*CD**2)/(ISA_density(h)*S*CL**3))
    
    Tr = Pr/V       #Convert power to thrust by dividing by airspeed
    return Tr   #returns the Power required 



"""Climb performance"""
#Varying thrust with velocity
def T_alt(T0,V,h):
    M = Mach(V,h)
    T = T0*0.369*M**-0.305
    return T
    
#Normal rate of climb at a given altitude and airspeed
def RC(Wcr,T0,V,S,A,e,CD0,h):
    k = 1. / (np.pi*A*e)
    T = T_alt(T0,V,h)
    RC = V*((T/Wcr) - 0.5*ISA_density(h)*V**2 * (Wcr/S)**(-1)*CD0 - (Wcr/S)*((2*k)/(ISA_density(h)*V**2)))
    return RC

#Maximum climb angle (or steepest climb) and corresponding airspeed and RC 
def steep_climb(T0,W,S,CD0,A,e,h):
    k = 1. / (np.pi*A*e)
    theta_max = np.arcsin((T0/W) - np.sqrt(4*CD0*(k))) 
    V_theta_max = np.sqrt((2/ISA_density(h))*(k/CD0)**0.5 * (W/S)*np.cos(theta_max))
    RC_max_theta = V_theta_max*np.sin(theta_max)
    
    return theta_max*(180/np.pi), V_theta_max,RC_max_theta   #return in degrees
    

#--------------------------------------MAIN PROGRAM------------------------------
"""Thrust required and available"""
V = np.arange(50,300,5)              #Criose velocity in m/s
H = np.arange(1000,12500,1000)       #Cruise altitude in m

min_Tr = []
min_Tr_M = []
plt.figure(1)  

for h in H:
    Tr_list = []
    M_list = []

    for v in V:
        Tr_list.append(Treq(Wcr,v,S,A,e,h,CD0)/1e03)
        M_list.append(Mach(v,h))
        
    plt.plot(M_list,Tr_list, label = "%s altitude" %h)
    
    k = Tr_list.index(min(Tr_list))
    min_Tr.append(min(Tr_list))
    min_Tr_M.append(M_list[k])
    
 
plt.plot(min_Tr_M,min_Tr,'ko',label = "Min. Treq" )  
plt.hlines(T0/1000,0,1.,"gray" ,'--', label = " Thrust available" )
plt.title("Power required" )
plt.xlabel("Mach number")
plt.ylabel("Required thrust [kN]")
plt.grid("True")
plt.legend()  

    
"""Climb rate"""
V = np.arange(50,300,5)              #Criose velocity in m/s
H = np.arange(1000,12500,1000)       #Cruise altitude in m

plt.figure(2)

RC_max = [] 
M_RC_max = []

for h in H:
    RC_list = []
    M_list = []

    for v in V:
        RC_list.append(RC(Wcr,T0,v,S,A,e,CD0,h))
        M_list.append(Mach(v,h))
        
    plt.plot(M_list,RC_list, label = "%s m" %h)
    

plt.title("Rate of climb" )
plt.xlabel("Mach number")
plt.ylabel("Rate of climb [m/s")
plt.legend()  
plt.grid("True")




    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




