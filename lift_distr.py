# -*- coding: utf-8 -*-
"""
Created on Sat May 11 13:30:31 2019

@author: floyd
"""
import numpy as np
import subprocess
import os
from matplotlib import pyplot as plt

def make_avl_file():
    # B777 used as reference aircraft
    S = 427.80
    span = 60.90
    MAC = 8.75
    AR = 8.67
    taper = 0.149
    qc_sweep = np.radians(31.60)
    dihedral = 0
    Cr = (2*S)/((1+taper)*span)
    print(Cr)
    Ct = Cr*taper
    chords = [Cr, Ct]
    CD_0 = 0.015
    Angle = 0.0
    
    dx = 0.25*Cr + span/2*np.tan(qc_sweep) - 0.25*Ct
    dz = span/2*np.tan(dihedral)
    
    x_loc_LE = [0, dx]
    y_loc_LE = [0, span/2]
    z_loc_LE = [0, dz]
    
    Ainc = [0.0, 0.0]
    spanwise_discretize_points = 50    #If you go too high then your computer is dead
    chordwise_discretize_point = 12     # " "
    
    with open("avl_testing.avl", "w") as text_file:
            print("Test Wing" +"\n"
            "#Mach" +"\n" + 
            str(0.7) +"\n"
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
                  "n2414.dat.txt", file=text_file)
make_avl_file()

def lift_distribution(CL):        
    p = subprocess.Popen(r"C:\Users\mathi\Documents\DSE\Bigger_is_Better\avl.exe", stdin=subprocess.PIPE, universal_newlines=True)
    set_CL = "a c " + str(CL)
    p.communicate(os.linesep.join(["load", "avl_testing","case", "mach0.7", "oper", set_CL, "x","fs", "endresult"]))          
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
    os.remove("endresult")
    return(elements)
output_avl = lift_distribution(0.8)

mac = 8.67
def get_correct_data(output_avl,mac):
    x_pos = []
    cl = []
    cd = []
    for i in range(len(output_avl)):
        x_pos.append(output_avl[i][1])
        cl.append(output_avl[i][4]/mac)
        cd.append(output_avl[i][8])
#    x_pos = x_pos[len(x_pos):int(len(x_pos)/2)-1:-1] + x_pos[0:int(len(x_pos)/2)]
#    cl = cl[len(x_pos):int(len(x_pos)/2)-1:-1] + cl[0:int(len(x_pos)/2)]
#    cd = cd[len(x_pos):int(len(x_pos)/2)-1:-1] + cd[0:int(len(x_pos)/2)]
#    plt.scatter(x_pos,cl)
#    plt.scatter(y_pos,cd)
#    plt.grid()
    return(x_pos,cl, cd)
x = get_correct_data(output_avl,mac)
print(x)