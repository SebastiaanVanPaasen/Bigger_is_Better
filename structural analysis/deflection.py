from loading_and_moment_diagrams import *
from math import *
from stress_distribution_wing import *

E = 69*10**9
I_yy = 0.026

G = 25.5*10**9
Am = 9.58926156*10**-2
t = 0.01
Fydistribution, Fzdistribution, Mydistribution, Mzdistribution, Tdistributionvalues, chord_section, HalfspanValues = load_diagrams(100)
    
airfoil = 'NACA3414.txt'
def ds_airfoil(airfoil):
    data = load_airfoil(airfoil)
    data_z = data[1]
    data_y = data[2]
    ds = 0
    for i in range(1,len(data_z)):
        ds_point = sqrt((data_z[i]-data_z[i-1])**2+(data_y[i]-data_y[i-1])**2)
        ds = ds + ds_point
    return(ds)

ds = ds_airfoil(airfoil)

def enclosed_properties(ds, chord_section, Am, t):
    ds_chord = list(np.array(chord_section)*ds)
    Am_chord = list(Am * (np.array(chord_section))*(np.array(chord_section)))
    t_chord = list(t * np.array(chord_section))
    shear_center_z = list(np.array(chord_section)*0.4)
    return ds_chord, Am_chord, t_chord, shear_center_z

ds_chord, Am_chord, t_chord, shear_center_z = enclosed_properties(ds, chord_section, Am, t)

def deflection_bending(Mzdistribution,HalfspanValues):
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

deflection_distribution = deflection_bending(Mzdistribution,HalfspanValues)

def angle_twist(Tdistributionvalues,HalfspanValues,airfoil,chord_section,Am,t):
    ds = ds_airfoil('NACA3414.txt')
    ds_chord, Am_chord, t_chord, shear_center_z = enclosed_properties(ds, chord_section, Am, t)
    angle_of_twist = [0.]
    angle_of_twist_per_cut = 0.
    for i in range(1,len(ds_chord)):
        angle_of_twist_section = ((Tdistributionvalues[i-1]*(HalfspanValues[i]-HalfspanValues[i-1]))/(4*(Am_chord[i]**2)*G))*(ds_chord[i]/t)
        angle_of_twist_per_cut = angle_of_twist_per_cut + angle_of_twist_section
        angle_of_twist.append(angle_of_twist_per_cut)
    return angle_of_twist, angle_of_twist_section


angle_of_twist, angle_of_twist_section = angle_twist(Tdistributionvalues,HalfspanValues,airfoil,chord_section,Am,t)

def deflection_torque(angle_of_twist,shear_center_z,chord_section, angle_of_twist_section):
    deflection_torque_LE = []
    deflection_torque_TE = []
    deflection_torque_LE_per_cut = 0.
    deflection_torque_TE_per_cut = 0.
    for i in range(len(angle_of_twist)):
        deflection_torque_LE_section = tan(angle_of_twist[i])*shear_center_z[i]
        deflection_torque_LE_per_cut = deflection_torque_LE_section + deflection_torque_LE_per_cut
        deflection_torque_TE_section = tan(angle_of_twist[i])*(shear_center_z[i]-chord_section[i])
        deflection_torque_TE_per_cut = deflection_torque_TE_section + deflection_torque_TE_per_cut
        deflection_torque_LE.append(deflection_torque_LE_section)
        deflection_torque_TE.append(deflection_torque_TE_section)
    return deflection_torque_LE, deflection_torque_TE


deflection_torque_LE, deflection_torque_TE = deflection_torque(angle_of_twist,shear_center_z,chord_section, angle_of_twist_section)

def deflection(deflection_torque_LE, deflection_torque_TE,deflection_bending):
    deflection_LE = list(np.array(deflection_torque_LE)+np.array(deflection_bending))
    deflection_TE = list(np.array(deflection_torque_TE)+np.array(deflection_bending))
    return deflection_LE, deflection_TE

deflection_LE, deflection_TE = deflection(deflection_torque_LE, deflection_torque_TE,deflection_distribution)

plt.plot(HalfspanValues,deflection_distribution)
#plt.plot(HalfspanValues,deflection_torque_TE)
plt.show()



#print(deflection_torsion)


#angle_of_twist = angle_twist(Tdistributionvalues,HalfspanValues,airfoil,chord_section,Am,t)

#deflection_bending = deflection(Mydistribution,HalfspanValues)
#deflection_torsion = deflection(Tdistributionvalues,HalfspanValues)

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