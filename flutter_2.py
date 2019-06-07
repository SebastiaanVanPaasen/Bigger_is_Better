from MOI import MOI, ell_lift, MOI2
import numpy as np
from matplotlib import pyplot as plt
b = 30 #SEMISPAN!!!!!!!
cl_cruise = 0.7
strutpos = 0#meter
wingtip_length = b - strutpos 
rootchord = 5
tipchord = 1
chordlength = 1#rootchord - strutpos/b*(rootchord-tipchord)  
airfoil_thickness = 0.001
sweep = np.radians(20) 
CL_alpha = 5#2*np.pi
S = 250

rho = .3
#V = 230
#q = .5*rho*V**2
#q = .3

#Material Properties
G = 26.9E9 #AL 7075-T6
E = 71.7E9 


#CROSS SECTIONAL PROPERTIES
momentsI = MOI2(chordlength, airfoil_thickness) 
I_theta = momentsI[3]
I_gamma = momentsI[0]
I_theta_gamma = momentsI[2]

K_theta = G*I_theta/wingtip_length #TORSIONAL STIFFNESS
K_gamma = E*I_gamma #BENDING STIFFNES
print("K_theta: ", K_theta)
print("K_gamma: ", K_gamma)

#Find lift distribution
lift_d, ypos = ell_lift(cl_cruise/2, b) #semi lift coefficient and semi span
index_ypos = np.where(ypos == strutpos)[0][0]
partial_lift = lift_d[index_ypos:]
mean = np.mean(partial_lift)
index_cp = np.argmin(abs(partial_lift-mean))
y_cp = ypos[index_ypos+index_cp]

e_cp = 0.01*chordlength
q = np.arange(0,100,.1)
values = []
V = np.sqrt(2*q/rho)
for x in q:
    
    a4 = I_theta * I_gamma - I_theta_gamma**2
    a2 = I_theta*(K_gamma+ x*S*CL_alpha*np.sin(sweep)*y_cp) + \
        I_gamma*(K_theta - x*S*CL_alpha*np.cos(sweep)*e_cp) - \
        I_theta_gamma*x*S*CL_alpha*(np.cos(sweep)*y_cp - np.sin(sweep)*e_cp)
    a0 = K_theta*K_gamma + K_theta*S*x*CL_alpha*np.sin(sweep)*y_cp - K_gamma*S*x*CL_alpha*np.cos(sweep)*e_cp 
    y = a2**2 - 4*a4*a0
    values.append(y)
#print(x)
#coef = [a4,a2,a0]
#roots = np.roots(coef)
plt.plot(V,values)
plt.grid()