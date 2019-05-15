from loading_and_moment_diagrams import *
from math import *


E = 1*10**9
I_yy = 5.6*(10**-1)

Fydistribution, Fzdistribution, Mydistribution, Mzdistribution, Tdistributionvalues, chord_section = load_diagrams(100)
    
airfoil = 'NACA3414.txt'
def ds(airfoil):
    data = load_airfoil(airfoil) 
    data_z = data[1]
    data_y = data[2]
    ds = 0
    for i in range(1,len(data_z)):
        ds_point = sqrt((data_z[i]-data_z[i-1])**2+(data_y[i]-data_y[i-1])**2)
        ds = ds + ds_point
    return(ds)
    
def enclosed_properties(ds, chord_section, Am, t):
    ds_chord = ds * chord_section
    Am_chord = Am * chord_section
    t_chord = t * chord_section
    return ds_chord, Am_chord, t_chord

def deflection(Mzdistribution,HalfspanValues):
    slope = 0.
    deflection = 0.
    slope_distribution = []
    deflection_distribution = []
    for i in range(len(HalfspanValues)-1):
        slope_distribution.append(slope)
        deflection_distribution.append(deflection)
        c1 = E*I_yy*slope
        c2 = E*I_yy*deflection
        slope = (Mzdistribution[i]*(HalfspanValues[i+1]-HalfspanValues[i])+c1)/(E*I_yy)
        deflection = ((Mzdistribution[i]*(HalfspanValues[i+1]-HalfspanValues[i])**2)+c1*(HalfspanValues[i+1]-HalfspanValues[i])+c2)/(E*I_yy)
    slope_distribution.append(slope)
    deflection_distribution.append(deflection)
    
    return deflection_distribution

def angle_twist(Tdistributionvalues,HalfspanValues,airfoil,chord_section,Am,t):
    ds = ds(airfoil)
    ds_chord, Am_chord, t_chord = enclosed_properties(ds, chord_section, Am, t)
    

deflection_bending = deflection(Mydistribution,HalfspanValues)
deflection_torsion = deflection(Tdistributionvalues,HalfspanValues)

#def enclosed_properties():
#
#
#
#
#def angle_of_twist(Tdistributionvalues,HalfspanValues,G,Am,ds,t):
#    return
#
#plt.plot(HalfspanValues,deflection_bending)
#plt.show()
#print(deflection_torsion)