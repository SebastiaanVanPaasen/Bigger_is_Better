# -*- coding: utf-8 -*-
"""
Created on Wed May  8 16:39:55 2019

@author: floyd
"""
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt


from planform import wing_parameters

M_cruise = 0.75      #inputs from different part 
surface_area = 200  #inputs from different part 
aspect_ratio = 9    #inputs from different part 
CL_cruise = 0.5
CD0 = 0.010
#variables in defining the boxwing geometry

delta_x_quarterchord = 30 #distance of quarter chords of upper and lower wing in case this increases
                            #then the upper wing is elongated as the lower wing is sized conventionally
height_boxwing_root = 10    #height of the boxed wing
angle_winglet = 60 #angle of the side winglets in degrees, increasing this angle then upper wing length decreases.

def boxed_wing_geometry(M_cruise, CL_cruise, surface_area, aspect_ratio, delta_x_quarterchord, height_boxwing_root, angle_winglet):
    #Use planform function to find span
    wing_geo = wing_parameters(M_cruise, CL_cruise, surface_area, aspect_ratio) #function from planform.py
    span = wing_geo[3]
    
    #Surface and aspect ratio single wing 
    surface_single = surface_area/2
    aspect_ratio_single = span**2/surface_single
    
    #FIND planform for lower wing
    wing_geo_lower = wing_parameters(M_cruise, CL_cruise, surface_single, aspect_ratio_single)
    sweep_lower = wing_geo_lower[0] #quarterchord sweep
    sweep_lower_LE = wing_geo_lower[1]
    
    #find height of the wingtip based on the height of the total boxwing and angle of the winglets.
    height_boxwing_tip = height_boxwing_root - span/2*np.tan(wing_geo_lower[6]) #ASSUMED THAT UPPER WING HAS NO DIHEDRAL ANGLE
    
    #longitudinal displacement of the winglet lower vs upper. 
    delta_x_winglet = height_boxwing_tip*np.tan(np.radians(angle_winglet))   #ASSUMED THAT WINGLET HAS NO TAPER
    
    #Upper sweep calc based on the delta x in quarterchord position
    sweep_upper = np.arctan(2*delta_x_quarterchord/span-np.tan(sweep_lower)-delta_x_winglet*2/span) 
    
    #taper calc based on upper sweep to obtain correct surface area
    taper_ratio_upper = 0.2*(2.-sweep_upper)
    tip_chord = wing_geo_lower[5]
    root_chord_upper = tip_chord/taper_ratio_upper
    
    x_tip_lower = span/2*np.tan(sweep_lower_LE)
    x_tip_upper = x_tip_lower + height_boxwing_tip*np.tan(np.radians(angle_winglet))
    x_root_upper = delta_x_quarterchord + 0.25*wing_geo_lower[4]-0.25*root_chord_upper
    z_tip_lower = span/2*np.tan(wing_geo_lower[6])
    
    x_loc_LE = np.array([0, x_tip_lower, x_tip_upper, x_root_upper])
    y_loc_LE = np.array([0, span/2, span/2,0])
    z_loc_LE = np.array([0,z_tip_lower, height_boxwing_root ,height_boxwing_root])
    chords = np.array([wing_geo_lower[4],tip_chord, tip_chord, root_chord_upper])
    #x_loc_TE = x_loc_LE + chords
    return(x_loc_LE,y_loc_LE,z_loc_LE,chords)

x_loc_LE, y_loc_LE, z_loc_LE, chords = boxed_wing_geometry(M_cruise, CL_cruise, surface_area, aspect_ratio, delta_x_quarterchord, height_boxwing_root, angle_winglet)

#OPTIONAL plotting of the leading edge
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.plot(x_loc_LE, y_loc_LE, z_loc_LE, color='black')
##ax.plot(x_loc_TE, y_loc_LE, z_loc_LE, color='black')
#plt.show()


def make_avl_file_bw(x_loc_LE, y_loc_LE, z_loc_LE, chords, M_cruise, surface_area, CD0, aspect_ratio):
    name = ["Lower wing", "Right winglet", "Upper wing"]
    Ainc = [0.0, 0.0, 0.0, 0.0] # delta to increase/decrease incidence angle locally
    Nspanwise = [0, 0, 0, 0]        #avl bullshit
    Sspace = [0, 0, 0, 0]           #same
    Angle = [0.0, 0.0, 0.0, 0.0] #Incidence angle
    
    
    
    with open("Output.avl", "w") as text_file:
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
        for j in range(3):
            print("\n" + "SURFACE" +"\n" +
                  str(name[j]), "\n"
                  "8  1.0  12  -2.0"+"\n"
                  "COMPONENT"+"\n" +
                  str(1), "\n"
                  "YDUPLICATE"+"\n" +
                  str(0.0), "\n" +
                  "ANGLE"+"\n"+
                  str(Angle[j]), file=text_file)
            for i in range(2):
                print("SECTION", file=text_file)            
                print(x_loc_LE[j+i],y_loc_LE[j+i],z_loc_LE[j+i],chords[j+i],Ainc[j+i],Nspanwise[j+i], Sspace[j+i], file=text_file)
            
make_avl_file_bw(x_loc_LE, y_loc_LE, z_loc_LE, chords, M_cruise, surface_area, CD0, aspect_ratio)
