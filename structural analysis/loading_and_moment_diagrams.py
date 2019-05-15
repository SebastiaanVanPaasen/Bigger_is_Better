# -*- coding: utf-8 -*-
"""
Created on Tue May 14 09:48:11 2019

@author: Mathilde
"""
import sys
sys.path.append("H:\DSE\Bigger_is_Better\Bigger_is_Better/class_I")
import numpy as np  ### Never use * to import stuff, as it makes it difficult to retrace where functions come from
import scipy as sp
import math as m
from scipy import interpolate  ### Useful to interpolate stuff
from scipy import integrate
from matplotlib import pyplot as plt
from lift_distr import *

##Moment Code


### Move the geometry definition to over here
CD0 = 0.02
S = 427.80  # m^2
AR = 8.67
taper = 0.149
Sweep0 = 0.558505  # rad
by = 20.  # m
b = 60.90
Cr = (S + np.tan(Sweep0) * by * (b / 4)) / (by + (b - by) * ((1 + taper) / 2))
Ct = Cr * taper
Cy = Cr - np.tan(Sweep0) * (by / 2)
Volume = 489.8545093  # m^3
Wing_W = 57461.507853787  # kg
Wing_Wf = 57461.507853787 * 9.81
specific_weight = Wing_Wf / Volume  # N/m^3
Fuel_W_tot = 1517632
Sweepsc = m.atan(m.tan(Sweep0) - 4 / AR * (0.4 * (1 - taper) / (1 + taper))) #sweep shear center
c_engine = 0 #position of start engine wrt chord
y_shear_center = 0.2   #position of y position shear center wrt
thrust_position = -1.5 # position of y of the thrust vector
x_fuel_begin = 0
x_fuel_end = 10.
start_eng_1 = 5.
start_eng_2 = 16.
n_engines = 4 
total_thrust = 1183376.56
engine_weight = 80067.989

# input for each critical case; as the lift distribution varies for each case
rho = 0.348331
V = 253
n = 2.5
W = 2000000.

def input_CL(n,S,V,rho,W):
    input_CL = 2.5*W/(0.5*rho*V**2*S)
    return input_CL
## Import File List:
output_avl = lift_distribution(input_CL(n,S,V,rho,W))
x_pos = get_correct_data(output_avl,c)[0]
x_total = x_pos[int((len(x_pos)/2)):][::-1] + x_pos[0:int((len(x_pos)/2))]

##Lift Code:
cl = get_correct_data(output_avl,c)[1]
cl_total = cl[int((len(x_pos)/2)):][::-1] + cl[0:int((len(x_pos)/2))]
PolyFitCurveCl = sp.interpolate.interp1d(x_total, cl_total, kind="cubic", fill_value="extrapolate")
# print('Lift Coefficients (highest order first, ending with 0th order term) are: \n{}\n'.format(PolyFitCurveCl))

##Drag Code:
cdi = get_correct_data(output_avl,c)[1]
cdi_total = cl[int((len(x_pos)/2)):][::-1] + cl[0:int((len(x_pos)/2))]
PolyFitCurveidrag = sp.interpolate.interp1d(x_total, cdi_total , kind='cubic', fill_value='extrapolate')


### Define your functions at the beginning of the program
def c(x):
    
    if x < (by / 2):
        c = Cr - 2 * x * ((Cr - Cy) / (by))
    if x > (by / 2):
        c = Cy - 2 * (x - (by / 2)) * ((Cy - Ct) / (b - by))
    return c

def S_cross_section(x):
    return c(x) * c(x) * 0.1

Vfuel = sp.integrate.quad(S_cross_section, x_fuel_begin, x_fuel_end)[0]
#print(Vfuel)
Fuel_W_tank = Fuel_W_tot/2
specific_W_f = Fuel_W_tank/Vfuel

def Loadcalculator(x0,Ff):
    Nloadcalculations = 100
    xnodevalues = np.linspace(x0, b / 2, Nloadcalculations)
    xleftvalues = xnodevalues[:-1]
    xrightvalues = xnodevalues[1:]
    xmiddlevalues = (xleftvalues + xrightvalues) / 2
    section_verticalforcelist = []
    section_horizontalforcelist = []
    Fz = 0
    Fy = 0
    Mz = 0
    My = 0
    T = 0
    W_f = 0
    D = 0
    L = 0
    Th = 0
    if x0 > start_eng_1:
        firstenginereachedyet = True
    else:
        firstenginereachedyet = False
    if x0 > start_eng_2:
        secondenginereachedyet = True
    else:
        secondenginereachedyet = False
        
    #print(xmiddlevalues)
    for i, x in enumerate(xmiddlevalues):
        ### Geometry calculations
        width = xrightvalues[i] - xleftvalues[i]
        chord = c(x)
        surfacearea = width * chord
        shear_center = 0.4 * c(x) + np.tan(Sweep0) * x  # assumption that the center of shear is at c/3 line down the wing span

        ###FINDING THE TOTAL FORCE DISTRIBUTION ALONG WING HALF SPAN

        ### Weight calculations
        t = 0.10
        section_volume = width * chord * t * chord
        section_weight = section_volume * specific_weight

        ### Lift calculations
        Cl = PolyFitCurveCl(x)
        section_lift = 0.5 * Cl * rho * (V ** 2) * surfacearea * n * 1.5* -1 # because it points in the (-)ive z-direction
        L = section_lift

        # print(surfacearea, Cl, section_lift)

        ###Fuel weight calculations
        # cruise_ff = 0.636572571
        # total_W_f = 0.5 * 1517632 * Ff * n * 1.5
        if x <= x_fuel_end:
            section_fuel_weight = specific_W_f * section_volume * Ff
        else:
            section_fuel_weight = 0

        ###Drag Calculations
        Cd = CD0 + PolyFitCurveidrag(x)
        section_drag = (0.5 * Cd * rho * (V ** 2) * surfacearea)* -1
#        print(section_drag)

        # (-) because it's in a direction opposite to thrust

        ###Thrust Calculations
        section_thrust = 0
        section_engineweight = 0
        if x > start_eng_1 and firstenginereachedyet == False and n_engines!=0:
            section_thrust = total_thrust / n_engines
            section_engineweight = engine_weight * n * 1.5
            Mz += -section_engineweight * (start_eng_1 - x0)
            My += section_thrust * (start_eng_1 - x0)
            firstenginereachedyet = True

            
        if x > start_eng_2 and secondenginereachedyet == False and n_engines!=0 and n_engines!=2:
            section_thrust = total_thrust / n_engines 
            section_engineweight = engine_weight * n * 1.5
            Mz += -section_engineweight * (start_eng_2 - x0)
            My += section_thrust * (start_eng_2 - x0)
            secondenginereachedyet = True
        
        ###Torque Calculations

        lift_position = (1 / 4) * c(x) + np.tan(Sweep0) * x  # assumtion that the lift along the span acts at c_0.25
        weight_position = (1 / 2) * c(x) + np.tan(Sweep0) * x  # assumption that the weight acts along the span acts at c_0.5
        fuel_position = (1 / 2) * c(x) + np.tan(Sweep0) * x  # assumption that the fuel acts along the span acts at c_0.5
        engine_position = c_engine * c(x) + np.tan(Sweep0) * x #position of the engine 

#        Cm = PolyFitCurveCm(x)
#        moment_aero = 0.5 * Cm * rho * (V ** 2) * (chord) * surfacearea * n * m.cos(Sweepsc)
        lift_torque = section_lift * (lift_position - shear_center) * m.cos(Sweepsc)
        weight_torque = section_weight * (weight_position - shear_center) * m.cos(Sweepsc)
        engine_torque = section_engineweight * (engine_position - shear_center) * m.cos(Sweepsc)
        fuel_torque = section_fuel_weight * (fuel_position - shear_center) * m.cos(Sweepsc)
        thrust_torque =  section_thrust * (y_shear_center - thrust_position)
        section_torque = lift_torque + weight_torque + engine_torque + fuel_torque + thrust_torque

        
        ###Net force sums
        section_verticalforceminusengine = section_lift + section_weight + section_fuel_weight
        section_verticalforce = section_verticalforceminusengine + section_engineweight
        section_horizontalforceminusthrust = section_drag
        section_horizontalforce = section_horizontalforceminusthrust + section_thrust

        ###Individual Force Distributions
        W_f += section_fuel_weight
        D += section_drag
        L += section_lift
        Th += section_thrust

        ###Net force calculations
        Fy += section_verticalforce
        Fz += section_horizontalforce

        ### Moment calculations
        Mz += -section_verticalforceminusengine * (x - x0)
        My += -section_horizontalforceminusthrust * (x - x0)

        ### Torque calculations
        T += section_torque
    
#    print("engine", section_engineweight)
#    print('L',L)
#    print('Cl',Cl)
    return Fy, Fz, Mz, My, L, W_f, D, Th, section_engineweight, T



def load_diagrams(N):  ### 100 nodes, so 99 beam elements
    
    HalfspanValues = np.linspace(0, b / 2 - 0.00001, N)
    Fydistribution = []
    Fzdistribution = []
    Mydistribution = []
    Mzdistribution = []
    Liftdistributionvalues = []
    Fuelweightdistributionvalues = []
    Dragdistributionvalues = []
    Thrustdistributionvalues = []
    Engine_distribution = []
    Tdistributionvalues = []
    chord_section = []
    
    for i, x in enumerate(HalfspanValues):
        Fy, Fz, Mz, My, L, W_f, D, Th, section_engineweight, T = Loadcalculator(x,1)
        chord = c(x)
        Fydistribution.append(Fy)
        Fzdistribution.append(Fz)
        Mzdistribution.append(Mz)
        Mydistribution.append(My)
        Liftdistributionvalues.append(L)
        Fuelweightdistributionvalues.append(W_f)
        Dragdistributionvalues.append(D)
        Thrustdistributionvalues.append(Th)
        Engine_distribution.append(section_engineweight)
        Tdistributionvalues.append(T)
        chord_section.append(chord)

   
#    HalfspanValues = np.linspace(0, b / 2 - 0.00001, N)
#    Fydistribution2 = []
#    Mydistribution2 = []
#    Mzdistribution2 = []
#    Liftdistributionvalues2 = []
#    Fuelweightdistributionvalues2 = []
#    Dragdistributionvalues2 = []
#    Thrustdistributionvalues2 = []
#    Engine_distribution2 = []
#    Fzdistribution2 = []
#    Tdistributionvalues2 = []
#    chord_section = []
#    
#    for i, x in enumerate(HalfspanValues):
#        Fy, Fz, Mz, My, L, W_f, D, Th, section_engineweight, T = Loadcalculator(x,0.6)
#        chord = c(x)
#        Fydistribution2.append(Fy)
#        Fzdistribution2.append(Fz)
#        Mzdistribution2.append(Mz)
#        Mydistribution2.append(My)
#        Liftdistributionvalues2.append(L)
#        Fuelweightdistributionvalues2.append(W_f)
#        Dragdistributionvalues2.append(D)
#        Thrustdistributionvalues2.append(Th)
#        Engine_distribution2.append(section_engineweight)
#        Tdistributionvalues2.append(T)
##        chord_section.append(chord)
#    
#    HalfspanValues = np.linspace(0, b / 2 - 0.00001, N)
#    Fydistribution3 = []
#    Mydistribution3 = []
#    Mzdistribution3 = []
#    Liftdistributionvalues3 = []
#    Fuelweightdistributionvalues3 = []
#    Dragdistributionvalues3 = []
#    Thrustdistributionvalues3 = []
#    Engine_distribution3 = []
#    Fzdistribution3 = []
#    Tdistributionvalues3 = []
#    chord_section = []
#    
#    for i, x in enumerate(HalfspanValues):
#        Fy, Fz, Mz, My, L, W_f, D, Th, section_engineweight, T = Loadcalculator(x,0)
##        chord = c(x)
#        Fydistribution3.append(Fy)
#        Fzdistribution3.append(Fz)
#        Mzdistribution3.append(Mz)
#        Mydistribution3.append(My)
#        Liftdistributionvalues3.append(L)
#        Fuelweightdistributionvalues3.append(W_f)
#        Dragdistributionvalues3.append(D)
#        Thrustdistributionvalues3.append(Th)
#        Engine_distribution3.append(section_engineweight)
#        Tdistributionvalues3.append(T)
##        chord_section.append(chord)
        
    plt.subplot(2,3,5)
    plt.subplot(2,3,1)
    plt.gca().set_title('Mz distribution')
    plt.plot(HalfspanValues, Mzdistribution)
    plt.subplot(2,3,2)
    plt.gca().set_title('Fz distribution')
    plt.plot(HalfspanValues, Fzdistribution)
    plt.subplot(2,3,3)
    plt.gca().set_title('My distribution')
    plt.plot(HalfspanValues, Mydistribution)
    plt.subplot(2,3,4)
    plt.gca().set_title('Fy distribution')
    plt.plot(HalfspanValues, Fydistribution)
    plt.subplot(2,3,5)
    plt.gca().set_title('T distribution')
    plt.plot(HalfspanValues, Tdistributionvalues)
    plt.show()

    return Fydistribution, Fzdistribution, Mydistribution, Mzdistribution, Tdistributionvalues, chord_section
Fydistribution, Fzdistribution, Mydistribution, Mzdistribution, Tdistributionvalues, chord_section = load_diagrams(100)
# TORQUE DISTRIBUTION CALCULATION
# M_max = max(Mydistribution)
# V_max = max (Fxdistribution)

# print(M_max, V_max)


#def DocumentAddition(Document,string1, string2):
#    Documentation = open(Document, 'a')
#    Documentation.write(string1+"   "+string2+"\n")
#    Documentation.close()
#
#    return
#""""
#DocumentAddition("Fz_data_2.5.txt", "Position Along Wing [m]", "Vertical Shear Force [N]")
#
#for i in range(len(HalfspanValues)):
#    DocumentAddition("Fz_data_2.5.txt",str(HalfspanValues[i]),str(Fzdistribution[i]))
#
#
#DocumentAddition("T_data_2.5.txt", "Position Along Wing [m]", "Torque [Nm]")
#
#for i in range(len(HalfspanValues)):
#    DocumentAddition("T_data_2.5.txt",str(HalfspanValues[i]),str(Tdistributionvalues[i]))
#
#
#DocumentAddition("Mx_data_2.5.txt", "Position Along Wing [m]", "Internal Moment [Nm]")
#
#for i in range(len(HalfspanValues)):
#    DocumentAddition("Mx_data_2.5.txt",str(HalfspanValues[i]),str(Mxdistribution[i]))
#
#
#"""
# plt.plot(HalfspanValues,Fxdistribution)
# plt.plot(HalfspanValues, Tdistributionvalues)

# fig = plt.figure()
# plt.plot(HalfspanValues,Mxdistribution)
# fig.suptitle('Internal Moment Distribution', fontsize=20)
# plt.xlabel('Position Along Wing Span', fontsize=18)
# plt.ylabel('Mx', fontsize=16)



#print(max(Mzdistribution))




# plt.plot(HalfspanValues, Engine_distribution)
# plt.plot(HalfspanValues,Liftdistributionvalues)
# plt.plot(HalfspanValues, Dragdistributionvalues)
# plt.plot(HalfspanValues, Thrustdistributionvalues)
# plt.plot(HalfspanValues,Fuelweightdistributionvalues)
#plt.figure(figsize=(8,8))
#plt.plot(HalfspanValues, Tdistributionvalues,label = 'Full fuel')
#plt.plot(HalfspanValues, Tdistributionvalues2,label = '60% fuel')
#plt.plot(HalfspanValues, Tdistributionvalues3,label = 'Zero fuel')
#legend = plt.legend(loc='upper center', shadow=True, fontsize = 12)
#plt.xlabel('$y_W [m]$', fontsize=12)
#plt.ylabel('$T [Nm]$', fontsize=12)
#
#plt.show()
#
#plt.close()
#plt.figure(figsize=(10,10))
#
#plt.plot(HalfspanValues, Mxdistribution,label = 'Full fuel')
#plt.plot(HalfspanValues, Mxdistribution2,label = '60% fuel')
#plt.plot(HalfspanValues, Mxdistribution3,label = 'Zero fuel')
#legend = plt.legend(loc='upper center', shadow=True, fontsize = 12)
#plt.xlabel('$y_W [m]$', fontsize=12)
#plt.ylabel('$M_x [Nm]$', fontsize=12)
#
#plt.show()
#
#plt.close()
#plt.figure(figsize=(10,10))
#
#plt.plot(HalfspanValues, Mzdistribution)
## plt.plot(HalfspanValues, Mzdistribution2,label = '60% fuel')
## plt.plot(HalfspanValues, Mzdistribution3,label = 'Zero fuel')
##legend = plt.legend(loc='upper center', shadow=True, fontsize = 12)
#plt.xlabel('$y_W [m]$', fontsize=12)
#plt.ylabel('$M_z [Nm]$', fontsize=12)
#
#plt.show()
#
#plt.close()
#plt.figure(figsize=(10,10))
#
#plt.plot(HalfspanValues, Fxdistribution)
#plt.xlabel('$y_W [m]$', fontsize=12)
#plt.ylabel('$F_x [N]$', fontsize=12)
#
#plt.show()
#
#plt.close()
#plt.figure(figsize=(10,10))
#
#plt.plot(HalfspanValues, Fzdistribution,label = 'Full fuel')
#plt.plot(HalfspanValues, Fzdistribution2,label = '60% fuel')
#plt.plot(HalfspanValues, Fzdistribution3,label = 'Zero fuel')
#legend = plt.legend(loc='upper center', shadow=True, fontsize = 12)
#plt.xlabel('$y_W [m]$', fontsize=12)
#plt.ylabel('$F_z [N]$', fontsize=12)
#
#plt.show()
#
#plt.close()

#print('in the order full fuel, 60%, emopty tank')
#
#print('Mx', max(Mxdistribution), max(Mxdistribution2), max(Mxdistribution3))
#print('Mz', max(Mzdistribution), max(Mzdistribution2), max(Mzdistribution3))
#print('Fz', max(Fzdistribution), max(Fzdistribution2), max(Fzdistribution3))
#print('Fx', max(Fxdistribution), max(Fxdistribution2), max(Fxdistribution3))
#print('T', max(Tdistributionvalues), max(Tdistributionvalues2), max(Tdistributionvalues3))


