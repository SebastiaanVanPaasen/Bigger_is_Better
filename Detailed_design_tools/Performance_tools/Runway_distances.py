# -*- coding: utf-8 -*-
"""
Created on Tue May 28 11:53:50 2019

@author: nikki
"""

#----------------------------IMPORT MODULES-----------------------------------
import numpy as np
import matplotlib.pyplot as plt


#------------------------------INPUTS------------------------------------------
"""B737-8 MAX inputs"""
MTOW = 82191*9.81
W_TO = MTOW
W_land = 69308*9.81

T_TO =  130*2*1000*1. 
#for land and take-off distanes times 0.7
#To account for normal thrust setting
T0 = T_TO

CL_maxto   =       2.9 
CL_max_land=       3.2
CD_land    =       0.09
CD_TO      =       0.07

A = 10.16
e = 0.7
S = 127

psi_TO     =       342.06     #Specific Thrust N/airflow [N/kg/s]
bypass     =       9       #Bypass ratio of the engine

g          =       9.80665
mu         =       0.03         #Ground roll friction on dry concrete/asphalt



"""Design input values"""
#MTOW       =      1520276.626
#W_TO       =       MTOW 
#MFW        =       170698.0674
#W_land     =       MTOW-0.8*MFW
#
#T_TO       =       398090.9781
#T0         =       2*216.17*1000 
#
#CL_maxto   =       2.9 
#CL_max_land=       3.2
#CD_land    =       0.09
#CD_TO      =       0.07
#
#A          =       15.       #Aspect ratio [-]
#e          =       0.6       #Oswald efficiency factor [-]
#S          =       211.1888478     #Surface area wing [m^2]
#
#psi_TO     =       342.06     #Specific Thrust N/airflow [N/kg/s]
#bypass     =       15       #Bypass ratio of the engine
#
#g          =       9.80665
#mu         =       0.03         #Ground roll friction on dry concrete/asphalt




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

    #Method 1 Torenbeek
    T_mean = 0.75*((5. + bypass)/(4.+ bypass))*T_TO  #quite low  
    mu_dash = 0.02 + 0.01*CL_maxto
    gamma_LOF = 0.9*(T_mean/W_TO) - (0.3/np.sqrt(A))
    h_to = 10.7                                      #[m] screen height of 35 ft
    
    S_run = (V_LOF**2/(2.*g))/((T_mean/W_TO)-mu_dash)
    S_air = (V_LOF**2)/(g*np.sqrt(2)) + (h_to/gamma_LOF)
    S_to = S_run + S_air
    
    #Method 2 Anderson
    kuc = 3.16e-05
    mu = 0.03
    Kt = (T_mean/W_TO) - mu
    G = ((16*h_to/34.3)**2) / ((1 + (16*h_to/34.3))**2) 
    k3 = (1/(np.pi*A*e))
    k1 = (1/3)*k3
    dCD0 = (W_TO/S)*kuc*(W_TO/9.81)**-0.215
    Ka = -(ISA_density(0)/(2*(W_TO/S)))*(0.02 + dCD0 + (k1 + (G/(np.pi*A*e))*CL_maxto**2) - mu*CL_maxto)
    
    sg = (1/(2*g*Ka))*np.log(1 + (Ka/Kt)*V_LOF**2) + 3*V_LOF
    R = (6.96*Vs**2)/9.81    
    theta = np.arccos(1- (h_to/R))
    sa = R*np.sin(theta)
    
    #Method 3: Flight mechanics course
    h_to = 10.7
    gamma = 3*np.pi/180  #climb angle
    
    Vs_TO = np.sqrt((W_TO*2)/(rho*S*CL_maxto))
    V_LOF = 1.05*Vs_TO
    CL_to = mu*np.pi*A*e
    
    T_mean = 0.75*((5. + bypass)/(4.+ bypass))*T_TO  #quite low      
    D_mean = 0.5*rho*(V_LOF/np.sqrt(2))**2*S*CD_TO
    L_mean =  0.5*rho*(V_LOF/np.sqrt(2))**2*S*CL_to
    
    Dg_mean = mu*(W_TO - L_mean)
    a_mean = (g/W_TO)*(T_mean - D_mean - Dg_mean)
    S_ground = V_LOF**2 / (2*a_mean)
    
    S_climb = (h_to - (1-np.cos(gamma))*(V_LOF**2/(0.15*g)))/ np.tan(gamma)
    S_trans = (V_LOF**2/(0.15*g))*np.sin(gamma)
    
    S_tot = S_ground + S_climb + S_trans
    
    return S_to,sa+sg,S_tot
    
 

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
    
    #Compute variables
    Vs = np.sqrt((W_TO*2.)/(S*rho*CL_maxto))    
    V2 = 1.2*Vs
    CL_to = mu*np.pi*A*e
    
    dgamma2 = (T_TO/W_TO)-((CL_to/CD_to)**(-1.)) - gamma2_min
    gamma_mean = 0.06 + (dgamma2)

    #Vx = V2*( ((1. + 2.*g*h_to/V2**2)/(1. + gamma_mean / (a_mean/g)))**0.5 - ((gamma_mean*g*(dt - 1.))/V2) )
    
    #Find overall equation
    S01 = Vx**2 / (2.*a_mean)                                    #Distance covered before engine failure at Vx
    S12 = (1./gamma_mean)*(((V2**2 - Vx**2)/(4.*a_mean)) + h_to) #Distance from engine failure up to save screen height at V2
    Sstop = (Vx**2/(2*a_stop)) + Vx*dt                           #If TO aborted (stop distance needed)


    S_continue = S01 + S12
    S_abord = S01 + Sstop

    return S_continue, S_abord   
    
    
def BFL(A,e,T_TO,W_TO,CD_to, CL_maxto, bypass,rho,g):      #Balanced field length: is when taking off is better than stopping in case of engine failure
    #Constants from torenbeek
    h_to = 10.7 #[m]
    gamma2_min = 0.024
    
    #Compute variables
    CL_to = mu*np.pi*A*e    
    dgamma2 = ((T_TO/W_TO)-(CL_to/CD_to)**(-1.)) - gamma2_min    
    CL2 = 0.694*CL_maxto
    T_mean = 0.75*((5. + bypass)/(4.+ bypass))*T_TO    
    mu_dash = 0.02 + 0.01*CL_maxto
    
    #Find overall equation
    a = 0.863/(1. + 2.3*dgamma2)    
    b = ((W_TO/S)/(rho*g*CL2)) + h_to    
    c = (1./(T_mean/W_TO - mu_dash)) + 4.6
    
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

def S_land(T0,g,W_land,S,rho,CL_max_land,CD_land):
    #Method 1 from torenbeek
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
    
    #Method 2 from Anderson
    mu = 0.03
    kuc = 3.16e-05
    Vf = 1.23*Vs
    
    R = Vf**2 / (0.2*g)
    hf = R*(1 - np.cos(3*np.pi/180.))
    theta = 2.*np.pi/180.
    sa = (15.3 - hf )/(np.tan(theta))
    sf = R*np.sin(theta)
    
    Vtd = 1.15*Vs
    k3 = (1/(np.pi*A*e))
    k1 = (1/3)*k3
    G = ((16*15.3/34.3)**2) / ((1 + (16*15.3/34.3))**2) 
    dCD0 = (W_land/S)*kuc*(W_land/9.81)**-0.215
    D = 0.5*ISA_density(0)*Vtd**2*CL_max_land*S
    Jt = (D/W_land) + mu
    Ja = (ISA_density(0)/(2*(W_TO/S)))*(0.02 + dCD0 + (k1 + (G/(np.pi*A*e))*CL_max_land**2) - mu*CL_max_land)
    
    sg = 1.5*Vtd + ((1/(2*9.81*Ja))*np.log(1 + (Ja/Jt)*Vtd**2))
    S_tor = sa+sf+sg

    #Method 3: Anderson
    Vs_land = np.sqrt((2*W_land)/(rho*S*CL_max_land))
    Vap = 1.3*Vs_land
    
    S_trans = 2.6*Vs_land
    
    gamma = 2.3*np.pi/180
    R = 1.3**2 * (((W_land/S)*(2/rho)*(1/CL_max_land))/(0.1*g))
    S_air = R*np.sin(gamma) + ((h_land - (1 - np.cos(gamma))*R)/np.tan(gamma))
    
    T_mean_rev = 0.25*((5. + bypass)/(4.+ bypass))*T0        
    D_mean = 0.5*rho*(Vap/np.sqrt(2))**2*S*CD_land
    L_mean =  0.5*rho*(Vap/np.sqrt(2))**2*S*CL_max_land*0.8
    
    S_brake = (W_land**2/(2*g*S)*(2/rho)*(1.3**2/CL_max_land**2)*(1. / (T_mean_rev + D_mean + mu*(W_land - L_mean))))
    
    return SL,S_tor, S_trans+S_air+S_brake



#------------------------------MAIN PROGRAM------------------------------------
rho = ISA_density(0)

#Standard take-off distance with all engines functioning
S_TO = TO_distance(W_TO,S,rho,CL_maxto,bypass,T_TO,A)  

#Take-off distance with engine failure: continued and abord
Vx = np.arange(0,76.,1)
S_TO_fail = TO_eng_fail(W_TO,g,S,rho,CL_maxto,A,e,T_TO,CD_TO,Vx)    
   
#Balenced field length
BFL = BFL(A,e,T0,W_TO,CD_TO, CL_maxto, bypass,rho,g)

S_land = S_land(T0,g,W_land,S,rho,CL_max_land,CD_land)

#Regulations according to flight mechanics
req_land = (10/6)*S_land[0]

if BFL > max(S_TO):
    req_TO = BFL
elif max(S_TO) > BFL:
    req_TO = max(S_TO)*1.15


#plt.hlines(S_land,Vx[0],Vx[-1],'k','--',label = "landing")
plt.hlines(BFL*1.01, Vx[0],Vx[-1],"gray","--",label = "BFL")
plt.plot(Vx,S_TO_fail[0],"g", label = "continued")
plt.plot(Vx,S_TO_fail[1],'r',label = "abord")
plt.xlabel("Engine failure speed [m/s]")
plt.ylabel("Distance covered [m]")
plt.title('Balenced field length (engine failure)')
plt.legend()

ax = plt.gca()
ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))

plt.show()    
    
print ("Standard TO length:" , S_TO, "m")
print ("Standard landing length:", S_land,"m")    
print ("Balanced field length:", BFL*1.01,"m")
print ()
print ("Required TO field length:", req_TO, "m")
print ("Required landing field length:", req_land, "m")
    
    
    
    
    





