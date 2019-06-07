import numpy as np
from matplotlib import pyplot as plt
from MOI import MOI2
#INPUTS

defl_without = 2375
defl_with = 477.8

m = 1000 #semi wingweight
airfoil_t = 0.001
rootchord = 5 
tipchord = 1
span = 30
chordlength = rootchord-0.75*(rootchord-tipchord)
I_theta = MOI2(chordlength,airfoil_t)
b = 0.127 #
e = 0.5
S_theta = 0.08587
K_h = 2818.8
K_theta = 37.3 
CL_alpha = 5#2*np.pi
S = 2*b

q_step = 1
q_start = 0
q_final = 150
q = np.arange(q_start,q_final,q_step)
nr_datapoints = int((q_final-q_start)/q_step)

a4 = np.full(nr_datapoints, m*I_theta- S_theta**2)
a2 = m*K_theta + I_theta*K_h - (2*m*e*b+S_theta)*q*CL_alpha
a0 = K_h*(K_theta - 2*e*b*q*S*CL_alpha)

coef = np.zeros((nr_datapoints,3))
sigma1 = []
sigma2 = []
sigma3 = []
sigma4 = []
omega1,omega2,omega3,omega4 = [],[],[],[]

for i in range(nr_datapoints):
    coef[i][0] = a4[i]
    coef[i][1] = a2[i]
    coef[i][2] = a0[i]
    x = [a4[i],0,a2[i],0,a0[i]]
    roots = np.roots(x).tolist()
    sigma1.append(np.real(roots[0]))
    sigma2.append(np.real(roots[1]))
    sigma3.append(np.real(roots[2]))
    sigma4.append(np.real(roots[3]))
    omega1.append(np.imag(roots[0]))
    omega2.append(np.imag(roots[1]))
    omega3.append(np.imag(roots[2]))
    omega4.append(np.imag(roots[3]))
plt.subplot(211)
plt.scatter(q,sigma1, color='C0')
plt.scatter(q,sigma2, color='C0')
plt.scatter(q,sigma3, color='C0')
plt.scatter(q,sigma4, color='C0')
plt.grid()

plt.subplot(212)
plt.scatter(q,omega1, color='C1')
plt.scatter(q,omega2, color='C1')
plt.scatter(q,omega3, color='C1')
plt.scatter(q,omega4, color='C1')
plt.grid()

## Mass Matrix
#M = np.array([[m, S_theta],[S_theta,I_theta]])
#
##Stiffness Matrix
#S = np.array([[K_h,0],[0,K_theta]])
#
##Aerodynamic Forces Matrix
#A = np.array([[0,-S*CL_alpha],[0,2*S*e*b*CL_alpha]])
