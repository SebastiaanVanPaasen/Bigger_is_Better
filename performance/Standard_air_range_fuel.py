# -*- coding: utf-8 -*-
"""
Created on Tue May 14 08:55:23 2019

@author: nikki
"""

#-------------------------------MODULES---------------------------------------
import numpy as np
import matplotlib.pyplot as plt

#-------------------------------DESIRED INPUTS-----------------------------------------
"""Weights"""
#Wcr        =       Aircraft weight in N during Cruise

"""Aicraft configuration"""
#CD0        =       Zero lift drag coefficient
#A          =       Aspect ratio
#e          =       Oswald efficiency factor 
#R_des      =       Design range [m]
#S          =       Surface area m^2

#------------------------STATISTICAL INPUTS----------------------------------
"""NOW THE SAME FOR REF. AND DESIGN AIRCRAFT, CHANGE WHEN ENGINE IS KNOWN"""
Ct0           = 12e-06      #Thrust Specific fuel conspumtion [kg/N/s] from B737MAX 
                            #proprtional to speed, divide by nominal conditions
                            #aka cruise speed
                            #12 future high bypass engines
                            #13.5 for bigger normal engines

#------------------------------ARBITRARY INPUTS-------------------------------

"""INPUTS: CHANGE ACCORDING TO DESIGN"""

MTOW = 82190.*9.81
OEW  = 45065*9.81
MLW  = 69308.*9.81
MZFW = 65952.* 9.81
MFW  = 20826.*9.81          # Maximum fuel weight (including reserve fuel)
W_fr = MFW/105. * 5.        #reserve fuel

A = 10.45
e = 0.85
CD0 = 0.020
g = 9.81
S = 152.4  #127 

             #[km]
Wcr = 0.8*MTOW              #63000*9.81#assumption for now
pax_max = 450
n = 1                       #load factor of number of passengers

#------------------------------VERIFICATION DATA--------------------------------

"""Create reference line of B737 8 Max"""
g = 9.81
MPW1 = 20882*g
MTOW1 = 82191*g
OEW1 = 45070*g
MFW1 = 31594*g


Wcr1 = 0.8*MTOW1             #ASSUMPTION!!! WHAT IS W_Cr ACTUALLY
#Wcr1 = 0.75*MTOW             #ASSUMPTION!!! WHAT IS W_Cr ACTUALLY

Mcr = 0.79
hcr = 12000                 #(m)
S1 = 124.5                    #m^2
b1 = 35.92                  #m
A1 = b1**2 / S1
e = 0.85
CD0 = 0.020

pax_ref = 200.
n_ref = 1.

#-----------------------------DEFINITIONS-------------------------------------
#Standard air range (SAR) = distances travelled per mass fuel burned
#Fuel burn is related to speed, altitude, thrust (or drag in steady flight)
#unit of SAR (m/s)/(kg/s)

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

 
          
def Mach(V,h):                  #enter V in km/h
    gamma = 1.4                 #enter h in m
    R = 287 #J/kg/K
    a = np.sqrt(gamma*R*ISA_temp(h))
    M = (V/3.6)/a
    return M 


def SAR(h,A,S,e,CD0,Ct0,Wcr):   #enter h in m, V in m/s
     V_list = np.linspace(600,1000,500)
     SAR = []
     for V in V_list:
        Ct = (Ct0/233.083)*(V/3.6)
        k = 1./(np.pi*A*e)
        q = 0.5*ISA_density(h)*(V/3.6)**2
        SARi = 1./((V/3.6) / ( (CD0 + k *(Wcr/(q*S))**2) *q*S*Ct ))     #in kg/m
        SAR.append(SARi*1000 )
        
     return SAR,V_list#in kg/km

   
#------------------------------MAIN PROGRAM------------------------------------

dh = 500                            #step size in altitude
H = range(7000,12500,dh)            #altitude range

min_SAR = []
V_minSAR = []

##For a given altitude (in def) run it for different speeds
#for h in H:   
#    SAR_list = SAR(h,A,S,e,CD0,Ct0,Wcr)[0]
#    V = SAR(h,A,S,e,CD0,Ct0,Wcr)[1]
#    print V
#    print SAR_list
#
#    min_SAR.append(min(SAR_list))
#    i = SAR_list.index(min(SAR_list))
#    V_minSAR.append(V[i])    
#    
#    for i in range(len(V)):         #Change velocity to Mach
#        V[i] = Mach(V[i],h)
#        
#    #plt.subplot(121)
#    plt.plot(V,SAR_list,label='%s altitude [m]' % h)
#    plt.title('Fuel consumption w.r.t. airspeed')
#    plt.xlabel("Mach number")
#    plt.ylabel("Fuel consumption [kg/km]") 
#
#
#For the reference case aim to stay below it:
SAR_ref = SAR(hcr,A1,S1,e,CD0,Ct0,Wcr1)[0]
V_ref = SAR(hcr,A1,S1,e,CD0,Ct0,Wcr1)[1]

#
#
for i in range(len(V_ref)):
    V_ref[i] = Mach(V_ref[i],hcr)   
    if 0.997*Mcr <= V_ref[i] <= 1.003* Mcr:
        SAR_ref_point = SAR_ref[i]


##PLot the single point of the ref. aircraft
##plt.plot(Mcr,SAR_ref_point,"mo", label = "Ref. aircraft")
##plt.legend()


#Plot with minimum SAR at different altitudes with corresponding Mach number
#for j in range(len(min_SAR)):
#    V_minSAR[j] = Mach(V_minSAR[j],H[j])
#        
#    plt.subplot(223)
#    plt.xlabel("Mach at minimum SAR ")
#    plt.ylabel("Minimum Fuel consumption [kg/km]")
#    plt.plot(V_minSAR[j],min_SAR[j],'o', label = 'altitude [m] %s ' % H[j])
#    plt.title('Minimum fuel consumption with corresponding Mach and altitude')
#    
##Plot reference aircraft    
#plt.plot(Mcr,SAR_ref_point,"mo", label = "Ref. aircraft")
#plt.legend(loc = "lower right")


#Fuel consumed per km per passenger/seat
for h in H:   
    pax = pax_max*n
    SAR_list = (SAR(h,A,S,e,CD0,Ct0,Wcr)[0])
    V = SAR(h,A,S,e,CD0,Ct0,Wcr)[1]

    for i in range(len(SAR_list)):
        SAR_list[i] = SAR_list[i]/pax_ref 
        
    for i in range(len(V)):                     #Change velocity to Mach
        V[i] = Mach(V[i],h)
    
    #plt.subplot(122)
    plt.plot(V,SAR_list,label='%s altitude [m]' % h)
    #plt.title('Fuel consumption per passenger w.r.t. Mach number')
    plt.xlabel("Mach number")
    plt.ylabel("Fuel consumption [kg/km/passenger]")

plt.plot(Mcr,SAR_ref_point/pax_ref,"mo", label = "Ref. aircraft")   
#plt.hlines(0.9*SAR_ref_point/pax_ref,0.5,1.,"gray",'--') 
plt.xlim(0.5,0.95)
plt.ylim(0.007,0.015)
plt.legend(loc = "upper left")

plt.show()

# print SAR_ref_point

#--------------------------------SENSITIVITY ANALYSIS-------------------------
#To change W_cr, A, CD0, Ct0 analyse a certain altitude
"""Inputs sensitivity analysis"""

#MTOW = 82190.*9.81
###OEW  = 45065*9.81
###MLW  = 69308.*9.81
###MZFW = 65952.* 9.81
###MFW  = 20826.*9.81          # Maximum fuel weight (including reserve fuel)
###W_fr = MFW/105. * 5.        #reserve fuel
###
#A =10.45
#Ct0           = 12e-06   
#e = 0.85
#CD0 = 0.020
#g = 9.81
#S = 127. 
#Wcr = 0.8*MTOW              #63000*9.81#assumption for now
###
#dh = 500                            #step size in altitude
#H = range(7000,12500,dh)            #altitude range


"""Altitude Analysis"""
#min_SAR = []
#V_minSAR = []
#
#"""Altitude Sensitivity analysis"""
#for h in H:   
#    SAR_list = SAR(h,A,S,e,CD0,Ct0,Wcr)[0]
#    V = SAR(h,A,S,e,CD0,Ct0,Wcr)[1]
#
#    min_SAR.append(min(SAR_list))
#    i = SAR_list.index(min(SAR_list))
#    V_minSAR.append(V[i])  
#    	
#diff_minSAR = []
#diff_V_minSAR = []  
# 
#for i in range(len(min_SAR)-1):
#    diff_minSAR.append(((min_SAR[i+1]-min_SAR[i])/min_SAR[i])*100.)
#    diff_V_minSAR.append(((V_minSAR[i+1]-V_minSAR[i])/V_minSAR[i])*100)
#    
#    
#print diff_minSAR
#print 
#print diff_V_minSAR
    


"""Variable Sensitivity analysis"""
#min_SAR_1 = []
#min_SAR_2 = []
#
#diff_SAR_tot = []
#diff_minSAR_tot = []
#diff_V_minSAR_tot = []
#
#for h in H:   
#    #calculate SAR for case 1 (base case)
#    SAR_list_1 = SAR(h,A,S,e,CD0,Ct0,Wcr)[0]
#    V_1 = SAR(h,A,S,e,CD0,Ct0,Wcr)[1]
#
#    min_SAR_1.append(min(SAR_list_1))
#    i = SAR_list_1.index(min(SAR_list_1))
#    V_minSAR_1 = (V[i])    
#    
#    for i in range(len(V)):         #Change velocity to Mach
#        V_1[i] = Mach(V_1[i],h)    
#    
#    #calculate SAR for case 2 (Change parameters)
#    SAR_list_2 = SAR(h,A,S,e,CD0,Ct0,Wcr)[0]
#    V_2 = SAR(h,A,S,e,CD0,Ct0,Wcr)[1]
#
#    min_SAR_2.append(min(SAR_list_2))
#    i = SAR_list_2.index(min(SAR_list_2))
#    V_minSAR_2=(V[i])    
#    
#    for i in range(len(V)):         #Change velocity to Mach
#        V_2[i] = Mach(V_2[i],h)
#        
#    #Compare the cases 
#    diff_SAR = []
#    for j in range(len(SAR_list_1)):
#        diff_SAR.append(((SAR_list_2[j]-SAR_list_1[j])/SAR_list_1[j])*100.)       #difference in percentage
#   #find average difference on SAR
#    diff_SAR1 = sum(diff_SAR)/len(diff_SAR)  
#    diff_SAR_min = (((min(SAR_list_2)-min(SAR_list_1))/min(SAR_list_1))*100.)
#    
#    diff_V_minSAR_tot.append(((V_minSAR_2-V_minSAR_1)/V_minSAR_1)*100)
#    diff_minSAR_tot.append(diff_SAR_min)
#    diff_SAR_tot.append(diff_SAR1)
#        
#print diff_minSAR_tot
#print
#print diff_V_minSAR_tot
#print
#print diff_SAR_tot
#        
    
    
#-------------------------------------------------------------------------------
#Other graph: Min SAR with respect to airspeed [km/h]
#for j in range(len(min_SAR)):
#    plt.subplot(222)
#    plt.xlabel("Airspeed at minimum SAR [km/h]")
#    plt.ylabel("Minimum Fuel consumption [kg/km]")
#    plt.plot(V_minSAR[j],min_SAR[j],'o', label = '%s altitude [m]' % H[j])
#    plt.title('Minimum fuel consumption with corresponding airspeed and altitude')
#  
#plt.legend()

"""Once speed and altitude are selected, more precies SAR can be made 
by taking into account the weight reduction due to fuel consumption"""





