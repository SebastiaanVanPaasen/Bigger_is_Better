import numpy  as np 
from matplotlib import pyplot as plt 
 
M_cruise = 0.75 
S = 238 
root_chord = 6.14 
span = 59.75 
CD_0 = 0.017697 
chordwise_discretize_point = 5 
spanwise_discretize_points = 8 
Angle = 0.0 
tip_chord = 1.826 
 
def define_radius(winglet_angle, winglet_radius,radius_points, winglet_LE_sweep): 
    xlim = np.cos(np.radians(winglet_angle)) 
    step = (xlim)/(radius_points) 
    y_point = (np.arange(0,xlim+step,step)) 
    z_point = -np.sqrt(1-y_point**2)+1 
    y_point,z_point = y_point*winglet_radius, z_point*winglet_radius 
#    print(x_point,y_point) 
#    plt.scatter(x_point,y_point)
    local_angle = np.degrees(np.arctan((y_point)/((radius+0.00000000001)-z_point))) 
    ratio = local_angle/90 
    location_along_radius = ratio*np.pi*radius/2 
    x_point = np.tan(np.radians(winglet_LE_sweep))*location_along_radius 
 
    return x_point,y_point,z_point, location_along_radius 
 
 
angle_w = 15 
radius = 1 
R_points = 4 
winglet_LE_sweep = 30 
 
x_loc_wingie,y_loc_wingie,z_loc_wingie, location_along_radius = define_radius(angle_w,radius,R_points, winglet_LE_sweep) 
 
#local_angle = np.degrees(np.arctan((y_loc_wingie)/((radius+0.00000000001)-z_loc_wingie))) 
#ratio = local_angle/90 
#location_along_radius = ratio*np.pi*radius/2 
#x_loc_wingie = np.tan(np.radians(winglet_LE_sweep))*location_along_radius 
 
A = 2.5 
taper_winglet = 0.3 
winglet_tip = tip_chord*taper_winglet 
winglet_b = A*.5*(tip_chord+winglet_tip) 
winglet_length = winglet_b-(np.pi*radius/2) 
 
 
tip_loc_y = y_loc_wingie[-1] + winglet_length*np.sin(np.radians(angle_w)) 
tip_loc_z = z_loc_wingie[-1] + winglet_length*np.cos(np.radians(angle_w)) 
tip_loc_x = winglet_b*np.tan(np.radians(winglet_LE_sweep)) 
wl_yloc = [y_loc_wingie[-1], tip_loc_y] 
wl_zloc = [z_loc_wingie[-1],tip_loc_z] 
wl_xloc = [x_loc_wingie[-1], tip_loc_x] 
 
#x_loc_wingie = np.zeros(np.size(y_loc_wingie)) 
#plt.scatter(x_loc_wingie,y_loc_wingie) 
#plt.axis('equal') 
 
 
chord_wingie = tip_chord - (location_along_radius/winglet_b)*(tip_chord-winglet_tip) 
wl_chord = [chord_wingie[-1], winglet_tip] 
 
with open("winglet.avl", "w") as text_file: 
    print("Winglet" + "\n" 
          "#Mach" + "\n" + 
          str(M_cruise) + "\n" 
          "#IYsym IZsym Zsym" + "\n" 
          "0  0  0" + "\n" 
          "#Sref  Cref  Bref" + "\n" + 
          str(S), str(round(root_chord, 3)), str(round(span)), "\n" 
          "#Xref Yref Zref" + "\n" 
          "0  0.0  0" + "\n" 
          "#CDcp" + "\n" + 
          str(round(CD_0, 3)), "\n" 
          "\n" + "SURFACE" + "\n" + 
          "Winglet", "\n" + 
          str(chordwise_discretize_point), "1.0 " + str(spanwise_discretize_points), "-2.0" + "\n" 
          "COMPONENT", "\n" + 
          str(1), "\n" 
#          "YDUPLICATE" + "\n" + 
#          str(0.0), "\n" + 
          "ANGLE" + "\n" + 
          str(Angle), file=text_file) 
    for i in range(np.size(x_loc_wingie)): 
        print("SECTION", file=text_file) 
        print(round(x_loc_wingie[i], 2), round(y_loc_wingie[i], 2), round(z_loc_wingie[i], 2), round(chord_wingie[i],2), str(0.0),str(0.0),str(0.0), file=text_file) 
           
    print("\n" + "SURFACE" + "\n" + 
      "Winglet2", "\n" + 
      str(chordwise_discretize_point), "1.0 " + str(spanwise_discretize_points), "-2.0" + "\n" 
      "COMPONENT", "\n" + 
      str(1), "\n" 
#          "YDUPLICATE" + "\n" + 
#          str(0.0), "\n" + 
      "ANGLE" + "\n" + 
      str(Angle), file=text_file) 
    for i in range(np.size(wl_xloc)): 
        print("SECTION", file=text_file) 
        print(round(wl_xloc[i], 2), round(wl_yloc[i], 2), round(wl_zloc[i], 2), round(wl_chord[i],2), str(0.0),str(0.0),str(0.0), file=text_file) 
 
#        print("CDCL" + "\n" 
#              "-1.54 0.002 0.69 0.0054 1.8 0.002", file=text_file) 
#    print("AFILE" + "\n" 
#                    "n2414.dat.txt", file=text_file)