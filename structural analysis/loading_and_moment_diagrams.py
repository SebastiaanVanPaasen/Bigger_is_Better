#sys.path.append("C:/Users/mathi/Documents/DSE/Bigger_is_Better")
#sys.path.append("C:/Users/Mels/Desktop/3e jaar TUDelft/DSE/code/Bigger_is_Better")
import numpy as np  ### Never use * to import stuff, as it makes it difficult to retrace where functions come from
import scipy as sp
import math as m

from matplotlib import pyplot as plt
from class_I.lift_distr import make_avl_file, lift_distribution, get_correct_data
import constants_and_conversions as cc
#from untitled1 import *
#import untitled1 as un
#from read_csv_input import read_output

#design = 'HIGH DD 2E STRUT'
#weights, wing, cruise_conditions = read_output(design)
#print(wing)
#print(weights)
#print(cruise_conditions)
### Move the geometry definition to over here
CD0 = 0.0177 
S = 238 
b = 60 
AR = 15 
taper = 0.297 
qcsweep= 0.5133
Sweep0 = m.atan(m.tan(qcsweep) - 4 / AR * (-0.25 * (1 - taper) / (1 + taper))) # rad

Cr = 6.14
Ct = Cr * taper
tc = 0.16
# Cy = 0.#Cr - np.tan(Sweep0) * (by / 2)

Wing_W = 205725 

Volume = un.V_tot  # m^3
specific_weight = Wing_W / Volume  # N/m^3
#Fuel_W_tot = 200716. 
density = un.Rho_fuel
Sweepsc = m.atan(m.tan(Sweep0) - 4 / AR * (0.4 * (1 - taper) / (1 + taper)))  # sweep shear center
c_engine = 0  # position of start engine wrt chord
thrust_position = -1.5  # position of y of the thrust vector
x_fuel_begin = 0.35 *(b / 2)
#x_fuel_end = 0.7 * (b / 2)
start_eng_1 = 7
start_eng_2 = 0.6 * (b / 2)
n_engines = 2
total_thrust = 386786  #396428.19#392632.62#418435.81 # 469612.93in Newton
engine_weight = 4100 * cc.g_0 + 25141 /n_engines #137279.1+25828.71#156890.4+25767.83#156890.4+21594.79#166696.05 + 23013.97 #137279.1+25828.71 in Newton
strut_position = b/2 - un.A_S_L[4]
strutforce = un.results[0][4]

# input for each critical case; as the lift distribution varies for each case
alt = 9000
rho = cc.Rho_0 * ((1 + (cc.a * alt) / cc.Temp_0) ** (-(cc.g_0 / (cc.R_gas * cc.a))))#0.32#0.32#0.43 #0.23
V = 235.42*1.3 #252.19#270.55#247.55#301.63
#print(V)
n_ult= 3.75
n = n_ult/1.5#4.7/1.5#4.21/1.5#4.4/1.5#4.4/1.5
W = 1536496 #- weights["W_F"]#1801946.31#1510125.47#1549762.26#1806203.58

def input_CL(S,V,rho,W):
    input_CL = W/(0.5*rho*V**2*S)
    return input_CL

print(input_CL(S, V, rho, W), "cl")

#print(input_CL(n, S, V, rho, W))
make_avl_file()
## Import File List:
output_avl = lift_distribution(input_CL(S, V, rho, W))
x_pos = get_correct_data(output_avl)[0]
#print(x_pos)
##Lift Code:
cl = get_correct_data(output_avl)[1]
PolyFitCurveCl = sp.interpolate.interp1d(x_pos, cl, kind="cubic", fill_value="extrapolate")

##Drag Code:
cdi = get_correct_data(output_avl)[2]
PolyFitCurveidrag = sp.interpolate.interp1d(x_pos, cdi, kind='cubic', fill_value='extrapolate')
#print(cdi)

### Define your functions at the beginning of the program
def c(z):
    c = Cr - ((Cr - Ct) / (b / 2)) * z
    return c

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
    W_sec = 0

    if x0 > start_eng_1:
        firstenginereachedyet = True
    else:
        firstenginereachedyet = False
    if x0 > start_eng_2:
        secondenginereachedyet = True
    else:
        secondenginereachedyet = False
    if x0 > strut_position:
        strutreachedyet = True
    else:
        strutreachedyet = False
##        
    #print(xmiddlevalues)
    for i, x in enumerate(xmiddlevalues):
        ### Geometry calculations
        width = xrightvalues[i] - xleftvalues[i]
        chord = c(x)
        surfacearea = width * chord
#        print("sur",surfacearea)
#        print('chord', chord)
#        print('width', width)
        shear_center = 0.4 * c(x) + np.tan(Sweep0) * x  # assumption that the center of shear is at c/3 line down the wing span

        ###FINDING THE TOTAL FORCE DISTRIBUTION ALONG WING HALF SPAN

        ### Weight calculations
        t = 0.14
        
        section_volume = un.Volumes[i]
        section_weight = section_volume * specific_weight * -1 

        ### Lift calculations
        Cl = PolyFitCurveCl(x)
        section_lift = 0.5 * Cl * rho * (V ** 2) * surfacearea * 1.5 * n  # because it points in the (-)ive z-direction
        
        ### Fuel calculations
        Fuel_W_tot = 200716./2 
        
        if x > x_fuel_begin and Fuel_W_tot > 0:
            section_fuel_weight = density * section_volume * cc.g_0 
        else:
            section_fuel_weight = 0

        Fuel_W_tot -= section_fuel_weight
            
        
#        if x <= x_fuel_end:
#            section_fuel_weight = specific_W_f * section_volume * Ff * -1
#        else:
#            section_fuel_weight = 0

        ###Drag Calculations
        Cd = CD0 + PolyFitCurveidrag(x) * n * 1.5
        section_drag = (0.5 * Cd * rho * (V ** 2) * surfacearea)*-1
    
        #print(section_drag)
        #        print(section_drag)

        # (-) because it's in a direction opposite to thrust
        section_strutforce = 0.
#        #Strut
        if x > strut_position and strutreachedyet == False:
            section_strutforce = strutforce *-1
            Mz += section_strutforce * (strut_position - x0)
            strutreachedyet = True

        ###Thrust Calculations
        section_thrust = 0
        section_engineweight = 0
        if x > start_eng_1 and firstenginereachedyet == False:
            section_thrust = total_thrust / n_engines
            section_engineweight = engine_weight * -1
 #           print(section_engineweight)
            Mz += section_engineweight * (start_eng_1 - x0)
            My += -section_thrust * (start_eng_1 - x0)
            firstenginereachedyet = True

        if x > start_eng_2 and secondenginereachedyet == False and n_engines != 0 and n_engines != 2:
            section_thrust = total_thrust / n_engines
            section_engineweight = engine_weight * -1
            Mz += section_engineweight * (start_eng_2 - x0)
            My += -section_thrust * (start_eng_2 - x0)
            secondenginereachedyet = True
        
        ###Torque Calculations

        lift_position = (1 / 4) * c(x) + np.tan(Sweep0) * x  # assumtion that the lift along the span acts at c_0.25
        weight_position = (1 / 2) * c(x) + np.tan(Sweep0) * x  # assumption that the weight acts along the span acts at c_0.5
        fuel_position = (1 / 2) * c(x) + np.tan(Sweep0) * x  # assumption that the fuel acts along the span acts at c_0.5
        engine_position = c_engine * c(x) + np.tan(Sweep0) * x  # position of the engine
        y_shear_center = 0.5*0.16*c(x)
        #        Cm = PolyFitCurveCm(x)
        #        moment_aero = 0.5 * Cm * rho * (V ** 2) * (chord) * surfacearea * n * m.cos(Sweepsc)
        strut_torque = -section_strutforce * (lift_position - shear_center) * m.cos(Sweepsc)
        lift_torque = -section_lift * (lift_position - shear_center) * m.cos(Sweepsc)
        weight_torque = -section_weight * (weight_position - shear_center) * m.cos(Sweepsc)
        engine_torque = -section_engineweight * (engine_position - shear_center) * m.cos(Sweepsc)
        fuel_torque = -section_fuel_weight * (fuel_position - shear_center) * m.cos(Sweepsc)
        thrust_torque = section_thrust * (y_shear_center - thrust_position)
        section_torque = lift_torque + weight_torque + engine_torque + fuel_torque + thrust_torque + strut_torque
        

        
        ###Net force sums
        section_verticalforceminusengine = section_lift + section_weight + section_fuel_weight
        section_verticalforce = section_verticalforceminusengine + section_engineweight
        section_horizontalforceminusthrust = section_drag
        section_horizontalforce = section_horizontalforceminusthrust + section_thrust

        ###Individual Force Distributions
        W_f += section_fuel_weight
        D += section_drag
        L += section_lift
#        print(L)
        Th += section_thrust
        
        W_sec += section_weight
        ###Net force calculations
        Fy += section_verticalforce + section_strutforce
        Fz += section_horizontalforce

        ### Moment calculations
        Mz += section_verticalforceminusengine * (x - x0)
        My += -section_horizontalforceminusthrust * (x - x0)

        ### Torque calculations
        T += section_torque
    

    return Fy, Fz, Mz, My, L, W_f, D, Th, section_engineweight, T, W_sec


#print(Loadcalculator(0, 1)[4])


def load_diagrams(N):  ### 100 nodes, so 99 beam elements

    HalfspanValues = np.linspace(0, b / 2 - 0.00001, N)
    Fydistribution = []
    Fzdistribution = []
    Mydistribution = []
    Mzdistribution = []
    Liftdistributionvalues = []
    Fuelweightdistributionvalues = []
    Sectionweightdistributionvalues = []
    Dragdistributionvalues = []
    Thrustdistributionvalues = []
    Engine_distribution = []
    Tdistributionvalues = []
    
    for i, x in enumerate(HalfspanValues):
        Fy, Fz, Mz, My, L, W_f, D, Th, section_engineweight, T, W_sec = Loadcalculator(x,1)
        Fydistribution.append(Fy)
        Fzdistribution.append(Fz)
        Mzdistribution.append(Mz)
        Mydistribution.append(My)
        Liftdistributionvalues.append(L)
        Fuelweightdistributionvalues.append(W_f)
        Sectionweightdistributionvalues.append(W_sec)
        Dragdistributionvalues.append(D)
        Thrustdistributionvalues.append(Th)
        Engine_distribution.append(section_engineweight)
        Tdistributionvalues.append(T)
    
    HalfspanValues = np.linspace(0, b / 2 - 0.00001, N)
    Fydistribution2 = []
    Mydistribution2 = []
    Mzdistribution2 = []
    Liftdistributionvalues2 = []
    Fuelweightdistributionvalues2 = []
    Sectionweightdistributionvalues2 = []
    Dragdistributionvalues2 = []
    Thrustdistributionvalues2 = []
    Engine_distribution2 = []
    Fzdistribution2 = []
    Tdistributionvalues2 = []
    
    for i, x in enumerate(HalfspanValues):
        Fy, Fz, Mz, My, L, W_f, D, Th, section_engineweight, T, W_sec = Loadcalculator(x,0.6)
        Fydistribution2.append(Fy)
        Fzdistribution2.append(Fz)
        Mzdistribution2.append(Mz)
        Mydistribution2.append(My)
        Liftdistributionvalues2.append(L)
        Fuelweightdistributionvalues2.append(W_f)
        Sectionweightdistributionvalues2.append(W_sec)
        Dragdistributionvalues2.append(D)
        Thrustdistributionvalues2.append(Th)
        Engine_distribution2.append(section_engineweight)
        Tdistributionvalues2.append(T)
    
    HalfspanValues = np.linspace(0, b / 2 - 0.00001, N)
    Fydistribution3 = []
    Mydistribution3 = []
    Mzdistribution3 = []
    Liftdistributionvalues3 = []
    Fuelweightdistributionvalues3 = []
    Sectionweightdistributionvalues3 = []
    Dragdistributionvalues3 = []
    Thrustdistributionvalues3 = []
    Engine_distribution3 = []
    Fzdistribution3 = []
    Tdistributionvalues3 = []
    
    for i, x in enumerate(HalfspanValues):
        Fy, Fz, Mz, My, L, W_f, D, Th, section_engineweight, T, W_sec = Loadcalculator(x,0)
        Fydistribution3.append(Fy)
        Fzdistribution3.append(Fz)
        Mzdistribution3.append(Mz)
        Mydistribution3.append(My)
        Liftdistributionvalues3.append(L)
        Fuelweightdistributionvalues3.append(W_f)
        Sectionweightdistributionvalues3.append(W_sec)
        Dragdistributionvalues3.append(D)
        Thrustdistributionvalues3.append(Th)
        Engine_distribution3.append(section_engineweight)
        Tdistributionvalues3.append(T)
    
#    print(Dragdistributionvalues)    
    plt.subplot(2,3,5)
    plt.subplot(2,3,1)
    plt.gca().set_title('Mz distribution')
    plt.plot(HalfspanValues, Mzdistribution)
    plt.ylim(-40000000,30000000)
    plt.xlabel('Position Along Wing Span [$m$]')
    plt.ylabel('Mz $[Nm]$')
    
    plt.subplot(2,3,2)
    plt.gca().set_title('Fz distribution')
    plt.plot(HalfspanValues, Fzdistribution)
    plt.ylim(-300000,300000)
    plt.xlabel('Position Along Wing Span [$m$]')
    plt.ylabel('Fz [$N$]')
    
    plt.subplot(2,3,3)
    plt.gca().set_title('My distribution')
    plt.plot(HalfspanValues, Mydistribution)
    plt.ylim(-3000000,500000)
    plt.xlabel('Position Along Wing Span [$m$]')
    plt.ylabel('My $[Nm]$')
    
    plt.subplot(2,3,4)
    plt.gca().set_title('Fy distribution')
    plt.plot(HalfspanValues, Fydistribution)
    plt.ylim(-1000000,5000000)
    plt.xlabel('Position Along Wing Span [$m$]')
    plt.ylabel('Fy [$N$]')
    
    plt.subplot(2,3,5)
    plt.gca().set_title('T distribution')
    plt.plot(HalfspanValues, Tdistributionvalues)
    plt.ylim(-1000000, 3000000)
    plt.xlabel('Position Along Wing Span [$m$]')
    plt.ylabel('T $[Nm]$')
    plt.show()
    
#    plt.plot(HalfspanValues,Liftdistributionvalues)
#    plt.show()
#    
#    plt.subplot(1,3,3)
#    plt.subplot(1,3,1)
##    plt.figure(figsize=(8,8))
##    plt.plot(HalfspanValues, Tdistributionvalues,dashes=[6, 2],label = 'Full fuel')
##    plt.plot(HalfspanValues, Tdistributionvalues2,dashes=[3, 1],label = '60% fuel')
#    plt.plot(HalfspanValues, Tdistributionvalues3,label = 'Zero fuel')
#    legend = plt.legend(loc='upper center', shadow=True, fontsize = 12)
#    plt.xlabel('Position Along Wing Span $[m]$', fontsize=12)
#    plt.ylabel('T $[Nm]$', fontsize=12)
#    
#    plt.subplot(1,3,2)
##    plt.figure(figsize=(10,10))
##    plt.plot(HalfspanValues, Mydistribution,dashes=[6, 2],label = 'Full fuel')
##    plt.plot(HalfspanValues, Mydistribution2,dashes=[3, 1],label = '60% fuel')
#    plt.plot(HalfspanValues, Mydistribution3,label = 'Zero fuel')
#    legend = plt.legend(loc='upper center', shadow=True, fontsize = 12)
#    plt.xlabel('Position Along Wing Span [$m$]', fontsize=12)
#    plt.ylabel('My [$Nm$]', fontsize=12)
#    
#    plt.subplot(1,3,3)
##    plt.figure(figsize=(10,10))
##    plt.plot(HalfspanValues, Mzdistribution, dashes=[6, 2], label = 'Full fuel')
##    plt.plot(HalfspanValues, Mzdistribution2,dashes=[3,1], label = '60% fuel')
#    plt.plot(HalfspanValues, Mzdistribution3, label = 'Zero fuel')
#    legend = plt.legend(loc='upper center', shadow=True, fontsize = 12)
#    plt.xlabel('Position Along Wing Span [$m$]', fontsize=12)
#    plt.ylabel('Mz [$Nm$]', fontsize=12)


#   plt.show()
#    print(Liftdistributionvalues[0])
    maxMz = max([max(Mzdistribution), max(Mzdistribution2), max(Mzdistribution3)])
    minMz = min([min(Mzdistribution), min(Mzdistribution2), min(Mzdistribution3)])
    maxMy = max([max(Mydistribution), max(Mydistribution2), max(Mydistribution3)])
    minMy = min([min(Mydistribution), min(Mydistribution2), min(Mydistribution3)])
    maxT = max([max(Tdistributionvalues), max(Tdistributionvalues2), max(Tdistributionvalues3)])
    maxFy = max([max(Fydistribution), max(Fydistribution2), max(Fydistribution3)])
    minFz = min([min(Fzdistribution), min(Fzdistribution2), min(Fzdistribution3)])
    maxFz = max([max(Fzdistribution), max(Fzdistribution2), max(Fzdistribution3)])

    return Mzdistribution, Mydistribution, Tdistributionvalues, maxMz,minMz, maxMy, minMy,maxT,maxFy,minFz, maxFz, Sectionweightdistributionvalues, Fuelweightdistributionvalues

#for j in range(len(strutposition)):
##    
#print(load_diagrams(100)[11])
#print(load_diagrams(100)[12])




