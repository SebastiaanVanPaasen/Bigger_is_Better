# plot stuff
import numpy as np
from conv_wing_avl import make_avl_file, find_clalpha
from matplotlib import pyplot as plt
M = [0.75,0.775,0.75,0.775,0.75]
#M = [0,0,0,0,0]
A = [13,13,14,15,13]
qc_sweep = np.array([0.453,0.513,0.453,0.513,0.453])
cr = np.array([6.15,6.04, 5.53,6.37, 5.85])
ct = np.array([1.9,1.8,1.71,1.89,1.81])
span = np.array([52.36,50.95,50.66,61.93,49.75]) 
taper = ct/cr

leading_edge_sweep = np.arctan(np.tan(qc_sweep) - (cr / (2 * span) * (taper - 1)))
half_chord_sweep = np.arctan(((ct / 2 + span / 2 * np.tan(leading_edge_sweep)) - (cr / 2)) / (span / 2))
S = [210.85, 199.70, 183.29, 255.70,190.37]
CD_0 = np.array([0.0234, 0.0242, 0.0259, 0.0211, 0.0239])
CD = np.array([0.0375, 0.0382, 0.0362, 0.0310, 0.0333])

cla = []
for i in range(len(M)):
    dihedral = np.deg2rad(3 - (np.rad2deg(qc_sweep[i]) / 10) - 2)

    make_avl_file(cr[i],ct[i],span[i],leading_edge_sweep[i],dihedral, S[i],CD_0[i],M[i],0,0,0,0,0,0,0,0,0,0,0,12, 5)
    x = (find_clalpha(M[i],CD_0[i])*57.3)
    cla.append(x)

#DATCOM CLALPHA
beta = np.sqrt(1-M[0]**2)
cla_datcom = 2*np.pi*A[0]/(2+np.sqrt(4+(A[0]*beta/0.95)**2*(1+np.tan(half_chord_sweep[0])**2/beta**2)))


#print(cla)
#print("HIGH 2E semi", "LOW 2E semidd", "HIGH 2E semidd strut", "HIGH DD STRUT", "High DD")
alpha = [-5,15]
alpha_0 = -3.171
for i in range(1):
    y0 = -alpha_0*cla[i]/57.3
    y1 = alpha[0]*cla[i]/57.3+y0
    y2 = alpha[1]*cla[i]/57.3+y0
    
    beta = np.sqrt(1-M[i]**2)
    cla_datcom = 2*np.pi*A[i]/(2+np.sqrt(4+(A[i]*beta/0.95)**2*(1+np.tan(half_chord_sweep[i])**2/beta**2)))
    print("Difference: " +str((cla_datcom-cla[i])/cla[i]*100))
    y1_datcom = alpha[0]*cla_datcom/57.3+y0
    y2_datcom = alpha[1]*cla_datcom/57.3+y0
    y = [y1,y2]
    y_datcom = [y1_datcom,y2_datcom]
    plt.plot(alpha,y_datcom)
    plt.plot(alpha,y)
plt.grid()
plt.legend(["AVL", "DATCOM", "LOW SDD", "HIGH SDD STRUT", "HIGH DD STRUT", "HIGH DD"])
plt.xlabel("Angle of Attack [deg]")
plt.ylabel("CL [-]")
plt.plot((-7,17),(0,0), color='black')
plt.xlim(-6,16)
plt.plot((0,0),(-.5,2.3), color='black')
plt.ylim(-0.5,2.3)
