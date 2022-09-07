from matplotlib import pyplot as plt
import numpy as np
x = []
y = []
chordlength = 1
file = open("C:/Users/floyd/OneDrive/Documenten/GitHub/Bigger_is_Better/avl/SC20612.txt",'r')
for line in file: 
    x.append(chordlength*(float(line.split()[0])))
    y.append(chordlength*(float(line.split()[1])))
#x_circle = np.arange(0,1.0,.01)
#y_circle = np.sqrt(1-x_circle**2)
y20 = np.asarray(y)/16*20
camber = []
camber20 = [] 
for i in range(int(len(x)/2-1)):
    camber.append((y[i]+y[-i])/2)
    camber20.append((y20[i]+y20[-i])/2)
    
#plt.plot(0.125*x_circle,0.125*y_circle)
#plt.plot(0.125*-x_circle, 0.125*y_circle)
plt.plot(x,y)
plt.plot(x,y20)
plt.plot(x[:101],camber, color='C0')
plt.plot(x[:101],camber20, color='C1')

#plt.axis('equal')
plt.ylim((-0.15,0.15))
plt.legend(['NACA SC(2)-0612','NACA SC(2)-0620' ])
plt.grid()