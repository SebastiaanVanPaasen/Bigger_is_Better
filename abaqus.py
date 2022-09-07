#from abaqus import *
#from abaqusConstants import * 
import math
#from matplotlib import pyplot as plt

local_chord1 = 5.928486308     #1.826
local_delta1 = 0  #15.02
local_chord2 = 1.854036494
local_delta2 = 14.627549886021889
strut_chord = 1
#m = mdb.models['Model-1'] 
#s = m.ConstrainedSketch(name='airfoil_sketch', sheetSize=15000.0) 
points = []
points2 = []
strut_airfoil = []
file = open("C:/Users/floyd/OneDrive/Documenten/GitHub/Bigger_is_Better/avl/SC20616.txt",'r')
file2 = open("C:/Users/floyd/OneDrive/Documenten/GitHub/Bigger_is_Better/avl/SC20612.txt",'r')
for line in file: 
    points.append([(float(line.split()[0])*local_chord1+local_delta1)*1000,(float(line.split()[1])*local_chord1)*1000])

for line2 in file2: 
    points2.append([(float(line2.split()[0])*local_chord2+local_delta2)*1000,(float(line2.split()[1])*local_chord2)*1000])

#    strut_airfoil.append([float(line.split()[0])*strut_chord*1000,float(line.split( )[1])*strut_chord*1000])
#s.Spline(points=points) 
#s.Spline(points=points2) 
#s.Spline(points=strut_airfoil)
for i in range(len(points)):
    plt.scatter(points[i][0],points[i][1], color='C0')
#plt.scatter((1826*0.4+15020),0, color='C1')
plt.scatter((local_chord1*0.4*1000),0,color='C1')
plt.axis('equal')
#plt.grid()
#
#height = 486.#mm
#x_c = 6140*0.4
#print(0.1*local_chord1*1000)