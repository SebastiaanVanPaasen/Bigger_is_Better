# -*- coding: utf-8 -*-
"""
Created on Tue May 14 10:31:45 2019

@author: floyd
"""
import numpy as np
import os
import subprocess
from run_conditions import define_run_condition
from matplotlib import pyplot as plt
from label_lines import *


#S = 32.5

#MAC = 8.75
#A = 10
#taper = 0.149
#qc_sweep = np.radians(35)
#dihedral = 0




def make_avl_file(S, A, taper, dihedral, qc_sweep):
    # B777 used as reference aircraft
    span = np.sqrt(S*A)
    Cr = (2*S)/((1+taper)*span)
    Ct = Cr*taper
    MAC = (2 / 3) * Ct * ((1 + taper + taper ** 2) / (1 + taper))
    chords = [Cr, Ct]
    CD_0 = 0.00
    Angle = 0.0
    
    dx = 0.25*Cr + span/2*np.tan(qc_sweep) + 0.25*Ct
    dz = span/2*np.tan(dihedral)
    
    x_loc_LE = [0, dx]
    y_loc_LE = [0, span/2]
    z_loc_LE = [0, dz]
    
    Ainc = [0.0, 0.0]
    spanwise_discretize_points = 12   #If you go too high then your computer is dead
    chordwise_discretize_point = 9    # " "
    
    with open("777.avl", "w") as text_file:
            print("Test Wing" +"\n"
            "#Mach" +"\n" + 
            str(0.0) +"\n"
            "#IYsym IZsym Zsym" +"\n"
            "0  0  0" +"\n" 
            "#Sref  Cref  Bref" +"\n"  +
            str(S), str(MAC), str(span), "\n"
            "#Xref Yref Zref" +"\n"
            "0  0.0  0" + "\n"
            "#CDcp" + "\n" + 
            str(round(CD_0,3)), "\n" 
            "\n" + "SURFACE" +"\n" 
            "Wing", "\n" + 
            str(chordwise_discretize_point), "1.0 " + str(spanwise_discretize_points), "-2.0"+"\n"
            "YDUPLICATE"+"\n" + 
            str(0.0), "\n" 
            "ANGLE"+"\n" +
            str(Angle), file=text_file)
            for i in range(2):
                print("SECTION", file=text_file)            
                print(round(x_loc_LE[i],3),round(y_loc_LE[i],3),round(z_loc_LE[i],3),round(chords[i],3),Ainc[i], file=text_file)        
            print("AFILE" + "\n"
                  "naca0010.dat", file=text_file)
    return(MAC)
#MAC = make_avl_file(S,A,taper,dihedral,qc_sweep)

def lift_distribution(CL,M_cruise, CD_0):   
   define_run_condition(M_cruise, CD_0)
   dir_path = os.path.dirname(os.path.realpath("lift_distr.py"))        
   p = subprocess.Popen(str(dir_path)+"/avl.exe", stdin=subprocess.PIPE, universal_newlines=True)
   set_alpha = "a c " + str(CL)
   p.communicate(os.linesep.join(["load", "777","case", "mach"+str(M_cruise)+".run", "oper", set_alpha, "x","fs", "endresult", "ft", "resulting_alpha"]))          
   lines = [line.rstrip('\n') for line in open('endresult')]
   elements = []
   count = 0
   for i in range(len(lines)):
       stripped = lines[i].split()
       if len(stripped) > 0 and stripped[0].isdigit():
           count = count + 1
           for j in range(len(stripped)):
               elements.append(float(stripped[j]))
   elements = np.reshape(np.array(elements),(count,-1))
   lines = [line.rstrip('\n') for line in open('resulting_alpha')]  
   CL = float(lines[23].split()[2])
   CDi = float(lines[25].split()[5])
   e = float(lines[27].split()[5])
   os.remove("endresult")
   os.remove("mach"+str(M_cruise)+".run")
   os.remove("resulting_alpha")
   return(elements, CL, CDi, e)
#CL = .4
#mach = .6
#cd0 = 0.000    
#output_avl = lift_distribution(alpha,mach,cd0)


def get_correct_data(output_avl,MAC,CD_0,A,j):
    y_pos = []
    cl = []
    cd = []
    for i in range(len(output_avl)):
        y_pos.append(output_avl[i][1])
        cl.append(output_avl[i][4]/MAC)
        cd.append(output_avl[i][8]+CD_0)
#    y_pos = y_pos[len(y_pos):int(len(y_pos)/2)-1:-1] #+ y_pos[0:int(len(y_pos)/2)]
    y_pos = y_pos[0:int(len(y_pos)/2)-1]
#    cl = cl[len(y_pos):int(len(y_pos)/2)-1:-1] #+ cl[0:int(len(y_pos)/2)]
    cl = cl[0:int(len(y_pos))]
#    cd = cd[len(y_pos):int(len(y_pos)/2)-1:-1] #+ cd[0:int(len(y_pos)/2)]
    plt.grid()
    plt.subplot(121)
    plt.plot(y_pos,cl, label= str(A[j]), color = "C0")
    plt.ylim((0,2.5))
    
#    plt.scatter(y_pos,cd)
#    plt.grid()
    return(y_pos,cl, cd)

#x = get_correct_data(output_avl,MAC,0.020)
"""------------------------------------------------------------------------------------------"""
 #ASPECT RATIO
#MAC = 8.75
#A = 10
#taper = 0.149
#qc_sweep = np.radians(35)
#dihedral = 0
#alpha = 3
#
#A = np.arange(5,21,2)
#CL_CDi = []
#CL_CDi_calc = []
#for j in range(len(A)):
#    MAC = make_avl_file(S,A[j],taper,dihedral,qc_sweep)
#    run_avl =lift_distribution(alpha,mach,cd0)
#    output_avl = run_avl[0]
#    CL_CDi.append(run_avl[1]/run_avl[2])
#    CL_CDi_calc.append(np.pi*A[j]*run_avl[3]/run_avl[1])
#    x = get_correct_data(output_avl,MAC,cd0,A,j)
#labelLines(plt.gca().get_lines(),zorder=2.5)    
#plt.subplot(122)
#plt.plot(A,CL_CDi, label= "avl induced drag")
##plt.figure(3)
#plt.plot(A,CL_CDi_calc, label= "approximated induced drag")
#plt.xlabel("Aspect Ratio")
#plt.ylabel("CL/CDi")
#plt.legend()
#plt.grid()
"""-------------------------------------------------------------------------------------"""
# #TAPER RATIO VS e
A = 5
S = 32.4
CL = 0.4739
mach = 0.01
cd0 = 0.000 

#MAC = 8.75
#A = 10
#taper = 0.149
#qc_sweep = np.radians(35)
dihedral = 0
taper = np.arange(0,1.05,0.05)
sweep = np.arange(0,35,5)
plt.figure(2)
plt.grid()
for k in range(len(sweep)):
    e = []
    for l in range(len(taper)):    
        MAC = make_avl_file(S,A,taper[l],dihedral,np.radians(sweep[k]))
        run_avl = lift_distribution(CL,mach,cd0)
        output_avl = run_avl[0]
        e.append(run_avl[3])
        
    print(e)
    #x = get_correct_data(output_avl,MAC,cd0,A,j)
    plt.plot(taper,e, label = str(sweep[k]))
labelLines(plt.gca().get_lines(),zorder=2.5)   
plt.show()


