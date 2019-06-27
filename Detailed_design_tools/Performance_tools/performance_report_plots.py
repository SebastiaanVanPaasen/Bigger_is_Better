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
Wto = 1497151.235
Wland = Wto- 0.9*170548.2285
Wcr = 1497151.235 - 0.4*170548.2285
T0 = 428.2*1000

Vcr = 218.7110308
S = 207.9561486
A = 13
e = 0.72
CD0 = 0.0233
CD = 0.0356#*2.77
CL_maxcr = 1.48
L_D = 17.43
L_Dmax= 18.32

T_land = T0*0.47
CD0_land = 0.0775

T_app = T0*0.397
T_TO = T0*0.385

CD0_TO = 0.0318


"""Needs new values for new design"""
CLmax_land = 3.05
CLmax_TO = 2.77
Vs_land = np.sqrt((2*Wland)/(1.225*S*CLmax_land))
Vs_TO = np.sqrt((2*Wto)/(1.225*S*CLmax_TO))

Va = 1.3*Vs_land
Vtd = 1.15*Vs_land
V_LOF = 1.2*Vs_TO 

print (1.2*V_LOF)

"""B737-8 MAX inputs"""
#Wto = 82191*9.81
#Wcr = (242670. - 0.4*20730)*9.81
#Wland = (242670. - 0.8*20730)*9.81
#T0 = 130*1000*2
#Vcr = 221.2
#S = 127
#A = 10.16
#e = 0.6
#CD0 = 0.017
#CD = 0.03
#
#T_land = T0*0.35
#T_app = T0*0.28
#T_TO = T0*0.37
#CLmax_land = 3.19
#CLmax_TO =2.89
#Vs_land = np.sqrt((2*Wland)/(1.225*S*CLmax_land))
#Vs_TO = np.sqrt((2*Wto)/(1.225*S*CLmax_TO))
#
#Va = 1.3*Vs_land
#Vtd = 1.15*Vs_land
#V_LOF = 1.2*Vs_TO 

#Upper and lower limits for the interpolation for the service ceiling
H_upper  = 12000
H_lower = 10000


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
    T = 0.9*T0*(ISA_density(h)/ISA_density(0))**m
    return T


"""Climb performance"""
#Normal rate of climb at a given altitude and airspeed (steady, no acceleration)
def RC(W,T0,V,S,A,e,CD0,h):
    k = 1. / (np.pi*A*e)
    T = T_alt(T0,h)
    RC = V*((T/W) - 0.5*ISA_density(h)*V**2 * (S/W)*CD0 - (W/S)*((2*k)/(ISA_density(h)*V**2)))
    return RC

def RC_unsteady(W,T,V,S,h,CD):  #at const. EAS thus accelerating during flight
    #T = T_alt(T0,h)
    M = Mach(V,h)
    vg_dvdh = 0.5668*M**2            #Constant EAS in tropospere
    D = 0.5*ISA_density(h)*V**2*S*CD
    C = ((T-D)*V)/(W*(1. + vg_dvdh))
    return C


#Maximum climb angle (or steepest climb) and corresponding airspeed and RC 
def steep_climb(T0,W,S,CD0,A,e,h):
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

"""Gliding unpowered descent"""
#V = np.arange(50,300,5)              #Criose velocity in m/s
#H = np.arange(1000,13000,1000)
#
#
#Vv_min_list = []
#range_list = []
#Vinf_list = []
#for i in range(len(H)):
#    Vv_min = -RD(Wcr,S,A,e,CD0,H[i])   
#    Vv_min_list.append(Vv_min)
#    V_inf = RD_V(Wcr,S,A,e,CD0,H[i])
#    Vinf_list.append(V_inf)
#    if H[i] == 9000:
#        print (glide_range(L_D,H[i])/1000)
#    range_list.append(glide_range(L_D,H[i])/1000)
#    
#    
#    
#     
#plt.figure(7)
#plt.plot(H,Vv_min_list)  
##plt.title("Minimum descend rate")
#plt.xlabel("Altitude [m]",fontsize='x-large'  )
#plt.ylabel("Rate of descent [m]",fontsize='x-large' )
#plt.grid(True)
#ax = plt.gca()
##ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
#ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
#
#plt.figure(8)
#plt.plot(H,range_list)  
#plt.vlines(9000,range_list[0],range_list[-1],"gray","--")
#plt.hlines(156.87,H[0],H[-1],"gray","--")
#plt.plot(9000,156.87,'ko')
##plt.title("Range during glide")
#plt.xlabel("Starting altitude [m]",fontsize='x-large' )
#plt.ylabel("Range [km]",fontsize='x-large' )
#plt.grid(True)
#ax = plt.gca()
##ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
#ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
#
#
#
#theta_min = (np.arctan(1/ (L_Dmax)))*(180/np.pi)
#print ("Minimum glide angle: ", theta_min," degrees")
#
#
#plt.show()









"""Steepest climb"""
#V = np.arange(50,300,5)              #Criose velocity in m/s
#H = np.arange(1000,13000,1000) 
#
#Vs_list = []
#theta_max_list = []
#V_thetamax_list = []
#RC_maxtheta_list = []
#Vmin_list = []
#
#for h in H: 
#    Vs = np.sqrt((2*Wcr)/(ISA_density(h)*S*CL_maxcr))             #to indentify whether V_theta_max is obtainable
#    Vs_list.append(Mach(Vs,h))                                  # Since often theta_max cannot be reached since req. V is lower than Vstall
#
#    theta_max = steep_climb(T0,Wcr,S,CD0,A,e,h)[0]
#    V_theta_max = steep_climb(T0,Wcr,S,CD0,A,e,h)[1]
#    RC_theta_max =  steep_climb(T0,Wcr,S,CD0,A,e,h)[2] 
#
#    theta_max_list.append(theta_max)
#    V_thetamax_list.append(Mach(V_theta_max,h))
#    RC_maxtheta_list.append(RC_theta_max)
#
#plt.figure(5)   
#
#
#plt.plot(H,theta_max_list)
##plt.title("Max. climb angle")
#plt.xlabel("Altitude [m]", fontsize = 'large')
#plt.ylabel("Climb angle [deg]", fontsize = 'large')
#plt.grid(True)
#ax = plt.gca()
##ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
#ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
#
#plt.figure(6)
#plt.plot(H,V_thetamax_list, label = "V at theta max")
#plt.plot(H,Vs_list,"--",label = "Stall speed")
##plt.title("Required speed at max. climb angle")
#plt.xlabel("Altitude [m]", fontsize = 'large')
#plt.ylabel("Mach number", fontsize = 'large')
#plt.grid(True)
#plt.legend(fontsize = 'large')
#
#ax = plt.gca()
##ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
#ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
#
#plt.figure(7)
#plt.plot(H,RC_maxtheta_list)
##plt.title("Rate of climb at max. climb angle")
#plt.xlabel("Altitude [m]", fontsize = 'large')
#plt.ylabel("Rate of climb [m/s]", fontsize = 'large')
#plt.grid(True)
#
#ax = plt.gca()
##ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
#ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))



"""Climb gradient"""  
#RC_app = RC(Wland,T_app,Va,S,A,e,CD0_land*0.95,0)
#RC_land = RC(Wland,T_land,Vtd,S,A,e,CD0_land,0)
#RC_to = RC(Wto,T_TO,V_LOF,S,A,e,CD0_TO,0) #at const. EAS thus accelerating during flight
#
#Vh_app = np.sqrt(Va**2 - RC_app**2)
#Vh_land = np.sqrt(Vtd**2 - RC_land**2)
#Vh_to = np.sqrt(V_LOF**2 - RC_to**2)
#
#theta_app = (RC_app/Vh_app*100)
#theta_land = RC_land/Vh_land*100
#theta_to = RC_to/Vh_to*100
#
#print (theta_app,theta_land,theta_to)  

        
"""Steady Climb rate"""
#Validation data 
#RC_validation = RC(Wto,T0,221.2,S,A,e,CD0,6096)
#RC_unsteady_vali =RC_unsteady(Wto,T0,221.2,S,6096,CD)
#print (RC_validation,RC_unsteady_vali)
#
#
#H = np.arange(1000,14000,1000)
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
#        RC_list.append(RC(Wto,T0,v,S,A,e,CD0,h))
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
#plt.legend(loc = "upper right",fontsize = "x-large" )  
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
#
#ax = plt.gca()
#ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
#ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
#
#plt.show()    
##  
#
"""Service and absolute ceilings"""
#Absolute ceiling is where RCmax = 0
#Service ceiling is where RDmax = 100 ft/min = 0.508 m/s (from Anderson book)
#H = list(H)
#k1 = H.index(H_lower)
#k2 = H.index(H_upper)
#
#serv_ceiling = H[k1] + ((0.508 - RC_max[k1])/(RC_max[k2]-RC_max[k1]))*(H[k2]-H[k1])
#abs_ceiling = H[k1] + ((0. - RC_max[k1])/(RC_max[k2]-RC_max[k1]))*(H[k2]-H[k1])
#
#print ("Service ceiling :",serv_ceiling,"m")
#print ("Absolute ceiling :",abs_ceiling,"m")
#
#plt.figure(3)
#
#plt.plot(RC_max,H)
#plt.vlines(0.508,H[0],H[-1],"gray","--",label = " Service ceiling")
#plt.hlines(serv_ceiling,RC_max[0],RC_max[-1],"gray","--")
#plt.vlines(0.,H[0],H[-1],"g","--",label = " Absolute ceiling")
#plt.hlines(abs_ceiling,RC_max[0],RC_max[-1],"g","--")
##plt.title("Service and absolute ceiling")
#plt.xlabel("Max. rate of climb [m/s]",fontsize='x-large')
#plt.ylabel("Altitude [m]",fontsize='x-large')
#plt.grid(True )
##plt.ylim(0,14000)
#plt.legend(fontsize='x-large')
#
#ax = plt.gca()
#ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
#ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
#

    
    
    
    
    
    
    
    
    




