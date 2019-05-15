# -*- coding: utf-8 -*-
"""
Created on Wed May  8 16:39:55 2019

@author: floyd
"""
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import os 
import subprocess
from run_conditions import define_run_condition

from planform import wing_parameters

M_cruise = 0.75      #inputs from different part 
surface_area = 200  #inputs from different part 
aspect_ratio = 10   #inputs from different part 
CL_cruise = 0.5
CD0 = 0.010
#variables in defining the boxwing geometry

delta_x_quarterchord = 30 #distance of quarter chords of upper and lower wing in case this increases
                            #then the upper wing is elongated as the lower wing is sized conventionally
height_boxwing_root = 6   #height of the boxed wing
angle_winglet = 45 #angle of the side winglets in degrees, increasing this angle then upper wing length decreases.
spanwise_discretize_points = 12
chordwise_discretize_points = 8
ROOT_dir = os.path.dirname(os.path.realpath("boxwing_geometry.py"))

def boxed_wing_geometry(M_cruise, CL_cruise, surface_area, aspect_ratio, delta_x_quarterchord, height_boxwing_root, angle_winglet):
    #Use planform function to find span
    wing_geo = wing_parameters(M_cruise, CL_cruise, surface_area, aspect_ratio, 0) #function from planform.py
    span = wing_geo[2]
    
    #Surface and aspect ratio single wing 
    surface_single = surface_area/2
    aspect_ratio_single = span**2/surface_single
    
    #FIND planform for lower wing
    wing_geo_lower = wing_parameters(M_cruise, CL_cruise, surface_single, aspect_ratio_single, 0)
    sweep_lower = wing_geo_lower[0] #quarterchord sweep
    sweep_lower_LE = wing_geo_lower[1]
    
    #find height of the wingtip based on the height of the total boxwing and angle of the winglets.
    height_boxwing_tip = height_boxwing_root - span/2*np.tan(wing_geo_lower[5]) #ASSUMED THAT UPPER WING HAS NO DIHEDRAL ANGLE
    
    #longitudinal displacement of the winglet lower vs upper. 
    delta_x_winglet = height_boxwing_tip*np.tan(np.radians(angle_winglet))   #ASSUMED THAT WINGLET HAS NO TAPER
    
    #Upper sweep calc based on the delta x in quarterchord position
    sweep_upper = np.arctan(2*delta_x_quarterchord/span-np.tan(sweep_lower)-delta_x_winglet*2/span) 
    
    #taper calc based on upper sweep to obtain correct surface area
    taper_ratio_upper = 0.2*(2.-sweep_upper)
    tip_chord = wing_geo_lower[4]
    root_chord_upper = tip_chord/taper_ratio_upper
    
    x_tip_lower = span/2*np.tan(sweep_lower_LE)
    x_tip_upper = x_tip_lower + height_boxwing_tip*np.tan(np.radians(angle_winglet))
    x_root_upper = delta_x_quarterchord + 0.25*wing_geo_lower[3]-0.25*root_chord_upper
    z_tip_lower = span/2*np.tan(wing_geo_lower[5])
    
    x_loc_LE = np.array([0, x_tip_lower, x_tip_upper, x_root_upper ])
    y_loc_LE = np.array([0, span/2, span/2, 0])
    z_loc_LE = np.array([0,z_tip_lower, height_boxwing_root ,height_boxwing_root])
    chords = np.array([wing_geo_lower[3],tip_chord, tip_chord, root_chord_upper ])
    #x_loc_TE = x_loc_LE + chords
    return(x_loc_LE,y_loc_LE,z_loc_LE,chords,sweep_upper)

x_loc_LE, y_loc_LE, z_loc_LE, chords,sweep_upper = boxed_wing_geometry(M_cruise, CL_cruise, surface_area, aspect_ratio, delta_x_quarterchord, height_boxwing_root, angle_winglet)

#OPTIONAL plotting of the leading edge
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.plot(x_loc_LE, y_loc_LE, z_loc_LE, color='black')
##ax.plot(x_loc_TE, y_loc_LE, z_loc_LE, color='black')
#plt.show()

def V_tail_sizing(root_chord_h, tip_chord_h, span_h, root_chord_v, tip_chord_v, span_v):
    S_v = .5*(root_chord_v+tip_chord_v)*span_v/2
    S_h = .5*(root_chord_h+tip_chord_h)*span_h
    Total_area = S_v + S_h
    angle = np.arctan(np.sqrt(S_v/S_h))
    return(Total_area, angle)

Total_area, angle = V_tail_sizing(3,1.5,13,2.5,1,11)

def V_tail_geometry(angle, Total_area, chords, y_loc_LE, sweep_upper, x_loc_LE):
    lateral_distance = height_boxwing_root * np.tan(angle) 
    lateral_ratio = lateral_distance/y_loc_LE[1]
    local_chord = chords[3] - (chords[3] - chords[2])*lateral_ratio
    forward_delta_x = lateral_distance*np.tan(sweep_upper)#-0.25*local_chord
    x_upper = x_loc_LE[3] - forward_delta_x
    
    pythagoras = np.sqrt(lateral_distance**2 + height_boxwing_root**2)
    x_lower = x_upper - pythagoras*np.tan(np.radians(46))
    root_chord = Total_area/pythagoras-local_chord 
    x_loc_LE_tail = [x_lower,x_upper]
    y_loc_LE_tail = [0,lateral_distance]
    z_loc_LE_tail = [0,height_boxwing_root]
    chords = [root_chord, local_chord]
    return(x_loc_LE_tail, y_loc_LE_tail,z_loc_LE_tail,chords)
geo_tail = V_tail_geometry(angle, Total_area, chords, y_loc_LE,sweep_upper,x_loc_LE)    

def make_avl_file_bw(x_loc_LE, y_loc_LE, z_loc_LE, chords, M_cruise, surface_area, CD0, aspect_ratio, geo_tail, spanwise_discretize_points, chordwise_discretize_points):
    name = ["Lower wing", "Right winglet", "Upper wing"]
    Ainc = [0.0, 0.0, 0.0, 0.0] # delta to increase/decrease incidence angle locally
    Nspanwise = [0, 0, 0, 0]        #avl bullshit
    Sspace = [0, 0, 0, 0]           #same
    Angle = [0.0, 0.0, 0.0, 0.0] #Incidence angle
        
    
    with open("prandtl.avl", "w") as text_file:
        print("Boxed Wing" +"\n"
              "#Mach" +"\n" + 
              str(M_cruise) +"\n"
              "#IYsym IZsym Zsym" +"\n"
              "0  0  0" +"\n" 
              "#Sref  Cref  Bref" +"\n" +
              str(surface_area), str(chords[0]), str(np.sqrt(surface_area*aspect_ratio)), "\n"          
              "#Xref Yref Zref" +"\n"
              "0.3  0.0  0.0806" + "\n"
              "CDcp" + "\n" +
              str(CD0), file=text_file)
        for j in range(2):
            print("\n" + "SURFACE" +"\n" +
                  str(name[j]), "\n" +
                  str(chordwise_discretize_points), "1.0 " + str(spanwise_discretize_points), "-2.0"+"\n"
                  "COMPONENT"+"\n" +
                  str(1), "\n"
                  "YDUPLICATE"+"\n" +
                  str(0.0), "\n" +
                  "ANGLE"+"\n"+
                  str(Angle[j]), file=text_file)
            for i in range(2):
                print("SECTION", file=text_file)            
                print(x_loc_LE[j+i],y_loc_LE[j+i],z_loc_LE[j+i],chords[j+i],Ainc[j+i],Nspanwise[j+i], Sspace[j+i], file=text_file)
        print("\n" + "SURFACE" +"\n" +
                  str(name[2]), "\n" + 
                  str(chordwise_discretize_points), "1.0 " + str(spanwise_discretize_points), "-2.0"+"\n"
                  "COMPONENT"+"\n" +
                  str(1), "\n"
                  "YDUPLICATE"+"\n" +
                  str(0.0), "\n" +
                  "ANGLE"+"\n"+
                  str(Angle[j]), file=text_file)   
        for i in range(2):
            print("SECTION", file=text_file)
            print(x_loc_LE[3-i],y_loc_LE[3-i],z_loc_LE[3-i],chords[3-i],Ainc[3-i],Nspanwise[3-i], Sspace[3-i], file=text_file)

        
        # VTAIL
        print("\n" + "SURFACE" +"\n" +
            "Stab", "\n" +
            str(chordwise_discretize_points), "1.0 " + str(spanwise_discretize_points), "-2.0"+"\n"
            "YDUPLICATE"+"\n" +
            str(0.0), "\n" +
            "ANGLE"+"\n"+
            str(0.0), file=text_file)
        for i in range(2):
            print("SECTION", file=text_file)            
            print(round(geo_tail[0][i],3),round(geo_tail[1][i],3),round(geo_tail[2][i],3),round(geo_tail[3][i],3),Ainc[i],Nspanwise[i], Sspace[i], file=text_file)        
            print("AFILE" + "\n"
              "n0010.dat.txt", file=text_file)
            
make_avl_file_bw(x_loc_LE, y_loc_LE, z_loc_LE, chords, M_cruise, surface_area, CD0, aspect_ratio, geo_tail, spanwise_discretize_points, chordwise_discretize_points)


def run_avl(cl_cruise,M_cruise, CD_0):
    define_run_condition(M_cruise, CD_0)        
    p = subprocess.Popen(str(ROOT_dir)  +"/avl.exe", stdin=subprocess.PIPE, universal_newlines=True)
    set_cl_cruise = "a c " + str(cl_cruise)
    p.communicate(os.linesep.join(["load", "prandtl","case", "mach"+str(M_cruise)+".run", "oper", set_cl_cruise, "x","ft", "endresult"]))          
    lines = [line.rstrip('\n') for line in open('endresult')]
    alpha = float(lines[15].split()[2])
    CD = float(lines[24].split()[2])
    e = float(lines[27].split()[5])
    print("Alpha is ", alpha)
    print("CD is ", CD)
    print("eff factor is ", e)
    os.remove("endresult")
    os.remove("mach"+str(M_cruise)+".run") 
run_avl(0.7,0.4,0.020)

def find_clalpha(M_Cruise,CD_0):
    define_run_condition(M_cruise, CD_0)  
    alpha_range = [0, 5]
    CL_range = []
    for j in range(len(alpha_range)):
        p = subprocess.Popen(str(ROOT_dir)  +"/avl.exe", stdin=subprocess.PIPE, universal_newlines=True)
        set_alpha = "a a " + str(alpha_range[j])
        p.communicate(os.linesep.join(["load", "prandtl","case", "mach"+str(M_cruise)+".run", "oper", set_alpha, "x","ft", "endresult"]))          
        lines = [line.rstrip('\n') for line in open('endresult')]
        CL = float(lines[23].split()[2])
        CL_range.append(CL)
        os.remove("endresult")
    CL_alpha = round((CL_range[1]-CL_range[0]) / (alpha_range[1] - alpha_range[0]),3)
    os.remove("mach"+str(M_cruise)+".run")
    return CL_alpha
print("CLa is ",find_clalpha(0.7,0.020))