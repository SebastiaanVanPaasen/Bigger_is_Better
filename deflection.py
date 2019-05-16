from loading_and_moment_diagrams import *

Fy, Fz, Mz, My, L, W_f, D, Th, section_engineweight, T = Loadcalculator(x,1)


    
N = 100  ### 100 nodes, so 99 beam elements
E = 1*10**9
I_yy = 5.6*(10**-1)
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

for i, x in enumerate(HalfspanValues):
    Fy, Fz, Mz, My, L, W_f, D, Th, section_engineweight, T = Loadcalculator(x,1)
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

deflection_bending = deflection(Mydistribution,HalfspanValues)
deflection_torsion = deflection(Tdistributionvalues,HalfspanValues)

def enclosed_properties():




def angle_of_twist(Tdistributionvalues,HalfspanValues,G,Am,ds,t):
    return

plt.plot(HalfspanValues,deflection_bending)
plt.show()
print(deflection_torsion)