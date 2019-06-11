import numpy as np
from matplotlib import pyplot as plt

K_theta = 1
e = 0.2
b = 1
E = 69E9 #MPa
S = 2*b
Cla = 2*np.pi
I_theta =5E-12#0.00443 #mm^4
m = 1
S_theta = m*e*S
K_h = E*I_theta
q = np.arange(0.01,1.01,0.01)
rho = .6
U = np.sqrt(2*q/rho)

a0 = K_h*(K_theta - 2*e*b*q*S*Cla)
a1 = q/U*S*Cla*K_theta
a2 = m*K_theta + I_theta*K_h - (2*m*e*b + S_theta)*q*S*Cla
a3 = q/U*S*Cla*(2*e*b*S_theta + I_theta)
a4 = m*I_theta - S_theta**2

p = []
q_scatter = []
for i in range(len(a0)):
    coef = [a4,a3[i],a2[i],a1[i],a0[i]]
    roots = np.roots(coef)
    for j in range(len(roots)):
        p.append(roots[j])
        q_scatter.append(q[i])
plt.scatter(q_scatter, np.real(p))
plt.grid()

#a = (a2**2-a4*a0)
#b = 1/(2*a4)+np.sqrt(-a2+a2**2-a4*a0)
#c = -a2/(2*a4)
#p1,p2,p3,p4 = [],[],[],[]
#
#for i in range(len(a)):
#    if a[i] > 0:
#        if b[i]>0:
#            p1.append( np.sqrt(b[i]) ) #pure real
#        else:
#            imag = np.sqrt(-b[i])
#            p2.append(complex(0,imag)) # pure imag
#    if a[i] < 0: 
#        if c[i]>0:
#            p3.append(complex(np.sqrt(c[i]),np.sqrt(-a[i]))) 
#        else:
#            p4.append(complex(0,np.sqrt(-c[i]+np.sqrt(-a[i]))))
#            
#plt.plot(q,p1)