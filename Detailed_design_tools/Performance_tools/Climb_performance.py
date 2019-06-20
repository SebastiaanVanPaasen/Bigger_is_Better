# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 15:36:13 2019

@author: nvanluijk
"""
import numpy as np
import matplotlib.pyplot as plt

#------------------------------INPUTS--------------------------------------------
MTOW       =       78220*9.81
MFW        =       31594*9.81
Wcr        =       MTOW-0.4*MFW 
T0         =       96.3*1000*2 #127.62*1000  #@ SEA LEVEL!!!
T_climb    =       127.9*1000*2

Vcr        =       233.

A          =       9.44       
e          =       0.85 
S          =       124.60 
CD0        =       0.020

CL_max     =       2.2
CL_maxcr   =       0.8
L_Dmax     =        20     
#CLcr       =       0.6
CDcr       =       0.05
#CD_TO      =       0.07
"""NOTE THIS VALUE DEPENDS ON THE ENGINE TYPE"""
m          =       1.3 #factor for the variation of thrust with alitude
 

h_climb = np.arange(1000,13000,1000)  #Altitudes to climb to from sea level
V_climb = np.arange(80.,225.,5)      #Velocities during the climb

#Upper and lower limits for the interpolation for the service ceiling
H_upper  = 12000
H_lower = 11000

#Altitude and velocity ranges for the plots
V = np.arange(50,300,5)              #Criose velocity in m/s
H = np.arange(5000,13000,1000)       #Cruise altitude in m

dh = 10000  #altitude that has to be descended

#Input needed for gliding range
#L/D fpr certain altitudes for now these values are assumed
L_D = [17,16.5,16,15.5,15,14.5,14,13.5,13,12.5,12,11.5,11]

#Wcr = 242670*9.81
#T0 = 341.5*2*1000   #kN
#S = 427.8
#A = 8.67
#e = 0.85
#CD0  = 0.026
#m = 1.3
#h = 3050.








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
    M = v/a
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

""" Max and min velocities""" 
def Vmax(Wcr,T0,S,A,e,CD0,h):
    #In order to minimise mistakes, split up eq. in more variables
    k = 1. / (np.pi*A*e)
    Tmax = T_alt(T0,h)
    Vmax = (((Tmax/Wcr)*(Wcr/S) + (Wcr/S)*np.sqrt((Tmax/Wcr)**2 - (4*CD0*k)))/(ISA_density(h)*CD0))**0.5
    return Vmax

def Vmin(Wcr,T0,S,A,e,CD0,h,CL_max,Vcr):
    Tmax = T_alt(T0,h)
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


"""Gliding flight: rate of descend"""
def RD(W,S,A,e,CD0,h):
    #Smalles rate of descend is obtained at max. (CL**(3/2))/CD
    k = 1. / (np.pi*A*e)
    CL_CD = 0.25*(3./(k*CD0**(1./3.)))**0.75
    Vv_min = - np.sqrt((2*W)/(ISA_density(h)*S))*(1./CL_CD)

    return Vv_min

def RD_V(W,S,A,e,CD0,h):
    k =  1. / (np.pi*A*e) 
    Vinf = np.sqrt((2./ISA_density(h))*np.sqrt((k*W)/(3*CD0*S)))
    return Vinf


def glide_range(L_D,dh):  
    R = (L_D)*dh
    return R




#--------------------------------------MAIN PROGRAM------------------------------
    


"""Thrust required and available"""
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
plt.grid(True)
plt.legend()  

"""Steady Climb rate"""
plt.figure(2)
plt.subplot(221)

RC_max = [] 
M_RC_max = []
V_RC_max = []
H_RC_max = []

for h in H:
    RC_list = []
    M_list = []
    V_list = []
    
    for v in V:
        RC_list.append(RC(Wcr,T0,v,S,A,e,CD0,h))
        M_list.append(Mach(v,h))
        V_list.append(v)
        
    plt.plot(M_list,RC_list, label = "%s m" %h)
    
    k = RC_list.index(max(RC_list))
    RC_max.append(max(RC_list))
    V_RC_max.append(V_list[k])
    M_RC_max.append(M_list[k])
    H_RC_max.append(h)

plt.plot(M_RC_max,RC_max," ko",label = " Max. RC" )
plt.title("Steady rate of climb" )
plt.xlabel("Mach number")
plt.ylabel("Rate of climb [m/s]")
plt.legend()  
plt.grid(True)


plt.subplot(222)
plt.plot(H,RC_max)
plt.title("Max. steady rate of climb" )
plt.xlabel("Altitude [m]" )
plt.ylabel("Rate of climb [m/s]" )
plt.grid(True)

plt.subplot(223)
plt.plot(H,M_RC_max)
plt.title("Mach number at Max. steady RC" )
plt.xlabel("Altitude [m]")
plt.ylabel("Mach number")
plt.grid(True)


"""Unsteady climb rate: acceleration or const. EAS"""
plt.figure(6)
plt.subplot(221)

RC_max_unst = [] 
M_RC_max_unst = []
V_RC_max_unst = []
H_RC_max_unst = []

for h in H:
    RC_list_unst = []
    M_list_unst = []

    for v in V:
        RC_list_unst.append(RC_unsteady(Wcr,T0,v,S,h,CDcr))
        M_list_unst.append(Mach(v,h))
        
    plt.plot(M_list_unst,RC_list_unst, label = "%s m" %h)
    
    k = RC_list_unst.index(max(RC_list_unst))
    RC_max_unst.append(max(RC_list_unst))
    V_RC_max_unst.append(v)
    M_RC_max_unst.append(M_list_unst[k])
    H_RC_max_unst.append(h)

plt.plot(M_RC_max_unst,RC_max_unst," ko",label = " Max. RC" )
plt.title("Unsteady rate of climb" )
plt.xlabel("Mach number")
plt.ylabel("Rate of climb [m/s]")
plt.legend()  
plt.grid(True)


plt.subplot(222)
plt.plot(H,RC_max_unst)
plt.title("Max. accelerated rate of climb" )
plt.xlabel("Altitude [m]" )
plt.ylabel("Rate of climb [m/s]" )
plt.grid(True)

plt.subplot(223)
plt.plot(H,M_RC_max_unst)
plt.title("Mach number at Max. unsteady RC" )
plt.xlabel("Altitude [m]")
plt.ylabel("Mach number")
plt.grid(True)

print ("RC max steady:", RC_max)
print ("RC max unsteady:", RC_max_unst)

"""Time to climb from sea level"""
#Normal time to climb 
plt.figure(4)


for h in h_climb:
    time = []
    climb_vel = []

    
    for v in V_climb:
        Rate = RC(MTOW,T_climb,v,S,A,e,CD0,h/2)
        time_climb = h/Rate

        time.append(time_climb/60)
        climb_vel.append(v)
 
    plt.plot(climb_vel,time,label = "%s To m" %h)
    
plt.title("Time to climb at steady RC" )
plt.xlabel("Climb velocity [m/s]")
plt.ylabel("Time to climb [min]")
plt.grid(True)
plt.legend()    


ax = plt.gca()
ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))

"""Service and absolute ceilings"""
#Absolute ceiling is where RCmax = 0
#Service ceiling is where RDmax = 100 ft/min = 0.508 m/s (from Anderson book)
H = list(H)
k1 = H.index(H_lower)
k2 = H.index(H_upper)

serv_ceiling = H[k1] + ((0.508 - RC_max[k1])/(RC_max[k2]-RC_max[k1]))*(H[k2]-H[k1])
abs_ceiling = H[k1] + ((0. - RC_max[k1])/(RC_max[k2]-RC_max[k1]))*(H[k2]-H[k1])

print ("Service ceiling :",serv_ceiling,"m")
print ("Absolute ceiling :",abs_ceiling,"m")

plt.figure(3)
plt.subplot(121)
plt.plot(RC_max,H)
plt.vlines(0.508,H[0],H[-1],"gray","--",label = " Service ceiling")
plt.hlines(serv_ceiling,RC_max[0],RC_max[-1],"gray","--")
plt.vlines(0.,H[0],H[-1],"g","--",label = " Absolute ceiling")
plt.hlines(abs_ceiling,RC_max[0],RC_max[-1],"g","--")
plt.title("Service and absolute ceiling")
plt.xlabel("Max. rate of climb [m/s]")
plt.ylabel("Altitude [m]")
plt.grid(True )
plt.legend()

ax = plt.gca()
ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))

Vs_list  = []
Vmax_list = []
for h in H:
    Vs = np.sqrt((2*Wcr)/(ISA_density(h)*S*CL_max))
    Vs_list.append(Vs)
    
    if h < serv_ceiling:
        Vmax_list.append(Vmax(Wcr,T0,S,A,e,CD0,h))
    else:
        Vmax_list.append(Vmax_list[-1])

    
plt.subplot(122)
plt.hlines(serv_ceiling,Vs_list[0],Vmax_list[0],"gray","--", label = "Service ceiling")
plt.plot(V_RC_max,H, label = "V @ RC max")
plt.plot(Vs_list,H, label = "Min. V limit")
plt.plot(Vmax_list,H,label = "Max. thrust limit")
plt.title("Flight envelope")
plt.grid(True)
plt.xlabel("Velocity [m/s]")
plt.ylabel("Altitude [m/s]")
plt.legend()    
ax = plt.gca()
ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))


"""Steepest climb"""
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


"""Gliding unpowered descent"""
Vv_min_list = []
range_list = []
Vinf_list = []
for i in range(len(H)):
    Vv_min = RD(Wcr,S,A,e,CD0,H[i])   
    Vv_min_list.append(Vv_min)
    V_inf = RD_V(Wcr,S,A,e,CD0,H[i])
    Vinf_list.append(V_inf)
    
    range_list.append(glide_range(L_D[i],H[i])/1000)
     
plt.figure(7)

plt.subplot(221)
plt.plot(H,Vv_min_list)  
plt.title("Minimum descend rate")
plt.xlabel("Altitude [m]" )
plt.ylabel("Rate of descent [m]" )
plt.grid(True)

plt.subplot(222)
plt.plot(H,range_list,' o' )  
plt.title("Range during glide")
plt.xlabel("Starting altitude [m]" )
plt.ylabel("Range [km]" )
plt.grid(True)

plt.subplot(223)
plt.plot(H,Vinf_list )  
plt.title("Velocity at minimum RD")
plt.xlabel("Altitude [m]" )
plt.ylabel("Velocity [m/s]" )
plt.grid(True)



ax = plt.gca()
ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))


theta_min = (np.arctan(1/ (L_Dmax)))*(180/np.pi)
print ("Minimum glide angle: ", theta_min," degrees")


plt.show()



#"""Hodograph"""
#plt.figure(6)
#
#for h in H:
#    Vv_list = []
#    Vh_list = []
#
#    for v in V:
#        Vv = RC(Wcr,T0,v,S,A,e,CD0,h)
#        Vh = np.sqrt(v**2  - Vv**2)
#  
#        Vv_list.append(Vv)
#        Vh_list.append(Vh)
#    
#    plt.plot(Vh_list,Vv_list, label = "%s m" %h)
#    
#plt.title("Hodograph climb performance")
#plt.xlabel("Horizontal velocity Vh [m/s]")
#plt.ylabel("Vertical velocity (RC) Vv [m/s]" )
#plt.grid(True)
#plt.legend()

#"""Min and max velocity depending on the Thrust available and drag"""
#Vmin_list = []
#Vmax_list = []
#Vs_list = []
#
#for h in H:
#    Vs = np.sqrt((2*Wcr)/(ISA_density(h)*S*CL_max))                         #Vstall with CLmax!!!
#    Vs_list.append(Mach(Vs,h))
#    
#    Vmin_list.append(Mach(Vmin(Wcr,T0,S,A,e,CD0,h,CL_maxcr,Vcr),h))         #Vmin during cruise!!!!!
#    Vmax_list.append(Mach(Vmax(Wcr,T0,S,A,e,CD0,h,Vcr),h))
#
#plt.figure(8)
#plt.plot(H,Vs_list,"--",label = "Stall speed with CL_max")
#plt.plot(H,Vmin_list, label = "Min. V in cruise")
#plt.plot(H,Vmax_list,label = "Max. V in cruise")
#plt.grid("True")
#plt.title("Min. d max. V depending on T and D")
#plt.xlabel("Altitude [m]")
#plt.ylabel("Mach number")
#plt.legend() 


    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




