# -*- coding: utf-8 -*-
"""
Created on Tue May 28 11:53:50 2019

@author: nikki
"""

#----------------------------IMPORT MODULES-----------------------------------
import numpy as np
import matplotlib.pyplot as plt


#------------------------------INPUTS------------------------------------------
"""Input values for the B737-800 aircraft"""
MTOW       =       78220*9.81 #Maximum take-off weight [N]
W_TO       =       MTOW       #Weight at take-off [N]
W_land     =       65310*9.81 #Maximum landing weight [N]

T_TO       =       0.2789*MTOW#Total static thrust of all engines at take-off [N]

CL_maxto   =       2.5  
CL_max_land=       3.3
CD_to      =       0.05 
CD_land    =       0.09

A          =       9.44       #Aspect ratio [-]
e          =       0.85       #Oswald efficiency factor [-]
S          =       124.60     #Surface area wing [m^2]

psi_TO     =       342.06     #Specific Thrust N/airflow [N/kg/s]
bypass     =       5.5        #Bypass ratio of the engine

g          =       9.80665

#------------------------------DEFINITIONS-------------------------------------
"""ISA definitions"""

def ISA_density(h):             # enter height in m
    M = 0.0289644               #kg/mol molar mass of Earth's air
    R = 8.3144590               #universal gas constant Nm/MolK
    
    if h < 11000:
        rho0 = 1.225            #kg/m^3
        T = 288.15              #K
        h0 = 0.                 #m
        a = -0.0065             #K/m
        rho = rho0 * (T/(T + a*(h-h0)))**(1. + ((g*M)/(R*a)))
        
    if h >= 11000:
        rho0 = 0.36391          #kg/m^3
        T = 216.65              #K
        h0 = 11000.             #m
        rho = rho0*np.e**((-g*M*(h-h0))/(R*T))
        
    return rho
    
    
def ISA_temp(h):
    if h < 11000:
        T = 288.15 - 0.0065*h   #in Kelvin
        return T
    if h >= 11000:
        return 216.65           #in Kelvin

def ISA_press(h):
    p0 = 101325.0
    T0 = 288.15 #[K]
    T = ISA_temp(h)
    a = -0.0065
    R = 8.3144590               #universal gas constant Nm/MolK
    
    if h <= 11000:
        p = p0*(T/T0)**(-g/(a*R))
    else:
        p=p0*np.e**(((-g/(R*T))*(h)))
        
    return p
 
          
def Mach(V,h):                  #enter V in m/s
    gamma = 1.4                 #enter h in m
    R = 287 #J/kg/K
    a = np.sqrt(gamma*R*ISA_temp(h))
    M = (V)/a
    return M 


"""Take-off distance: all engines functioning"""
def TO_distance(W_TO,S,rho,CL_maxto,bypass,T_TO,A): #Required distance to pass screen height (30 ft) at a speed of 1.3*VV_stall
    Vs = np.sqrt((W_TO*2.)/(S*rho*1.3*CL_maxto))   
    V_LOF = 1.2*Vs
    
    T_mean = 0.75*((5. + bypass)/(4.+ bypass))*T_TO  #quite low  
    mu_dash = 0.02 + 0.01*CL_maxto
    gamma_LOF = 0.9*(T_mean/W_TO) - (0.3/np.sqrt(A))
    h_to = 10.7                                      #[m] screen height of 35 ft
    
    S_run = (V_LOF**2/(2.*g))/((T_mean/W_TO)-mu_dash)
    S_air = (V_LOF**2)/(g*np.sqrt(2)) + (h_to/gamma_LOF)
    S_to = S_run + S_air
    
    return S_to
    
 

"""TO with engine failure"""
#gamma2_min =       Correction factor
#           =       0.024; 0.027; 0.030
#           =       2; 3; 4 number of engines

#dt         =       4.5 [s] reaction time
#mu         =       Ground friction coefficient
#           =       0.02 for concrete

def TO_eng_fail(W_TO,g,S,rho,CL_maxto,A,e,T_TO,CD_to,Vx):
    #Constants from torenbeek
    a_stop = 0.37*g
    a_mean = 0.25*g
    h_to = 10.7 #[m]
    gamma2_min = 0.024
    dt = 4.5
    mu = 0.02
    
    #Conpute variables
    Vs = np.sqrt((W_TO*2.)/(S*rho*CL_maxto))    
    V2 = 1.2*Vs
    CL_to = mu*np.pi*A*e
    
    dgamma2 = (T_TO/W_TO)-((CL_to/CD_to)**(-1.)) - gamma2_min
    gamma_mean = 0.06 + (dgamma2)

    #Vx = V2*( ((1. + 2.*g*h_to/V2**2)/(1. + gamma_mean / (a_mean/g)))**0.5 - ((gamma_mean*g*(dt - 1.))/V2) )
    
    #Find overall equation
    S01 = Vx**2 / (2.*a_mean)                                 #Distance covered before engine failure at Vx
    S12 = (1./gamma_mean)*(((V2**2 - Vx**2)/(2.*g)) + h_to)   #Distance from engine failure up to save screen height at V2
    Sstop = (Vx**2/(2*a_stop)) + Vx*dt                        #If TO aborted (stop distance needed)


    S_continue = S01 + S12
    S_abord = S01 + Sstop
    
    return S_continue, S_abord   
    
    
def BFL(A,e,T_TO,W_TO,CD_to, CL_maxto, bypass,rho,g):      #Balanced field length: is when taking off is better than stopping in case of engine failure
    #Constants from torenbeek
    h_to = 10.7 #[m]
    gamma2_min = 0.024
    mu = 0.02
    
    #Compute variables
    CL_to = mu*np.pi*A*e    
    dgamma2 = ((T_TO/W_TO)-(CL_to/CD_to)**(-1.)) - gamma2_min    
    CL2 = 0.694*CL_maxto
    T_mean = 0.75*((5. + bypass)/(4.+ bypass))*T_TO    
    mu_dash = 0.02 + 0.01*CL_maxto
    
    #Find overall equation
    a = 0.863/(1. + 2.3*dgamma2)    
    b = ((W_TO/S)/(rho*g*CL2)) + h_to    
    c = (1./(T_mean/W_TO - mu_dash)) + 2.7
    
    BFL = a*b*c
    return BFL
    
   
"""Landing distance"""
#Inputs
#gamma_mean =       During landing is approx. 0.1
#h_land     =       15.3 [m] screen height for landing
#dn         =       0.1 (change in load factor during landing)
#a_stop/g   =       0.35-0.45 without thrust reversers
#           =       0.40 - 0.50 thrustreversers and spoilers
#           =       0.50 - 0.60 + nose wheel braking

def S_land(g,W_land,S,rho,CL_max_land):
    #Constants from torenbeek
    h_land = 15.3
    dn = 0.1
    gamma_mean = 0.08
    a_stop = 0.55*g
    
    #Compute variables   
    Vs = np.sqrt((W_land*2.)/(S*rho*CL_max_land)) 
    Va = 1.4*Vs
    Vtd = np.sqrt(Va**2*(1. - (gamma_mean**2/dn)))    

    #Find equation
    S_air = (1./gamma_mean)*(((Va**2 - Vtd**2)/(2.*g)) + h_land)
    S_run = Vtd**2 / (a_stop)

    SL = S_air + S_run
    return SL    

#------------------------------MAIN PROGRAM------------------------------------
rho = ISA_density(0)

#Standard take-off distance with all engines functioning
S_TO = TO_distance(W_TO,S,rho,CL_maxto,bypass,T_TO,A)  

#Take-off distance with engine failure: continued and abord
Vx = np.arange(0,76.,1)
S_TO_fail = TO_eng_fail(W_TO,g,S,rho,CL_maxto,A,e,T_TO,CD_to,Vx)    
   
#Balenced field length
BFL = BFL(A,e,T_TO,W_TO,CD_to, CL_maxto, bypass,rho,g)


S_land = S_land(g,W_TO,S,rho,CL_max_land)

plt.hlines(S_land,Vx[0],Vx[-1],'k','--',label = "landing")
plt.hlines(BFL, Vx[0],Vx[-1],"gray","--",label = "BFL")
plt.plot(Vx,S_TO_fail[0],"g", label = "continued")
plt.plot(Vx,S_TO_fail[1],'r',label = "abord")
plt.xlabel("Engine failure speed [m/s]")
plt.ylabel("Distance covered [m]")
plt.title('Balenced field length (engine failure)')
plt.legend()
plt.show()    
    
print (" Standard TO length:" , S_TO, "m")
print (" Standard landing length:", S_land,"m")    
    
    
    
    
    
    
    





