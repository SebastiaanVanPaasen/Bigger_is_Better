import numpy  as np 
import os
import subprocess
from run_conditions import define_run_condition
from matplotlib import pyplot as plt
# GENERAL INPUTS
M_cruise = 0.72 
S = 207.9561486
root_chord = 5.89
span = 52
CD_0 = 0.0233-0.01
#CD_0 = CD_0 - 0.00243
chordwise_discretize_point = 8 
spanwise_discretize_points = 12 
Angle = 0.0 
tip_chord = 2.1
#span = 59.75302768
taper =0.357
qc_sweep = 0
LE_sweep = np.arctan(np.tan(qc_sweep) - (root_chord / (2 * span) * (taper- 1)))
dihedral = np.radians(-1.5)
intermed_chord = root_chord - 0.25*(root_chord-tip_chord)

#WINGLET INPUTS
#angle_w = 15 
radius = .5 
R_points = 3 
winglet_LE_sweep = 30 
A = 2.5
taper_winglet = 0.3 

#ENGINE INPUTS:
engine_xloc = 16.867
engine_yloc = 7.0
engine_zloc = -2.0
engine_placement = [engine_xloc,engine_yloc,engine_zloc]

#STRUT
strutpos = 18 #meter
wing_loc = 19.085
local_chord = root_chord - 20/(span/2)*(root_chord-tip_chord)
x_strut = [strutpos*np.tan(LE_sweep)+wing_loc+0.25*local_chord, strutpos*np.tan(LE_sweep)+wing_loc+0.25*local_chord]
y_strut = [0,strutpos]
oz_strut = [-3.34*2, (span/2*np.tan(dihedral))-.8]
strut_chords = [.99,.99]

def engine_avl(engine_placement):
    engine_d = 3.15
    engine_l = 5.60
    y_engine = np.array([0.0,0.5,0.866,1.0,0.866,0.5,0.0,-0.5,-0.866,-1,-0.866,-0.5,0.0])*engine_d/2 + engine_placement[1]
    z_engine = np.array([1.0,0.866,0.5,0,-0.5,-0.866,-1.0,-0.866,-0.5,0.0,0.5,0.866,1.0])*engine_d/2 + engine_placement[2]
    engine_chords = np.full(np.size(y_engine),engine_l)
    x_engine = np.full(np.size(y_engine), engine_placement[0])
    return(x_engine, y_engine, z_engine, engine_chords)


def avl(root_chord, tip_chord,span, LE_sweep, dihedral, angle_w,radius,R_points, winglet_LE_sweep, taper_winglet, A, intermed_chord, x_engine, y_engine, z_engine, engine_chords):
    wing_loc = 19.08509676
    #DEFINITIONS FOR MAIN WING
    chords = [root_chord, intermed_chord ,tip_chord]
    Ainc = [0.0,0.0,0.0]
    Nspanwise = [0, 0, 0]
    Sspace = [0, 0, 0]
    x_loc_LE = [0 + wing_loc,0.25*span / 2 * np.tan(LE_sweep) + wing_loc,span / 2 * np.tan(LE_sweep)+ wing_loc]
    y_loc_LE = [0, 0.25*span/2,span / 2]
#    print(x_loc_LE[2]-x_loc_LE[0])
    z_loc_LE = [0, 0.25*span/2*np.tan(dihedral), span / 2 * np.tan(dihedral)] 
#    print(x_loc_LE[2] - x_loc_LE[0])
    span_v  = 6.08
    h_chords = [4.73,1.89]
    x_loc_h = 50.37
    span_h = 13.24
    h_sweep = np.radians(15.57845313)
    x_loc_LE_h = [x_loc_h,  span_h / 2 * np.tan(h_sweep)+ x_loc_h]
    y_loc_LE_h = [0, span_h/2 ]
    z_loc_LE_h = [span_v, span_v]

    
    x_loc_v = 44.38
    v_sweep = np.radians(44.52928643)
    x_loc_LE_v =[x_loc_v, span_v  * np.tan(v_sweep)+ x_loc_v]
#    print(x_loc_LE_v)
    y_loc_LE_v = [0,0]
    z_loc_LE_v = [0, span_v]
    v_chords = [4.92, 3.20]

    
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
     
    
     
    x_loc_wingie,y_loc_wingie,z_loc_wingie, location_along_radius = define_radius(angle_w,radius,R_points, winglet_LE_sweep) 
    
    x_loc_wingie = x_loc_wingie + x_loc_LE[-1] 
    y_loc_wingie = y_loc_wingie + y_loc_LE[-1]
    z_loc_wingie = z_loc_wingie + z_loc_LE[-1]

    def winglet_tip(tip_chord,taper_winglet, A,x_loc_wingie,z_loc_wingie,y_loc_wingie,angle_w,winglet_LE_sweep, location_along_radius):
        winglet_tip = tip_chord*taper_winglet 
#        print(winglet_tip)
        winglet_b = A*.5*(tip_chord+winglet_tip) 
        winglet_length = winglet_b-(np.pi*radius/2) 
        
#        print(winglet_length)         
         
        tip_loc_y = y_loc_wingie[-1] + winglet_length*np.sin(np.radians(angle_w)) 
        tip_loc_z = z_loc_wingie[-1] + winglet_length*np.cos(np.radians(angle_w)) 
        tip_loc_x = x_loc_wingie[-1] + winglet_b*np.tan(np.radians(winglet_LE_sweep)) 
        wl_yloc = [y_loc_wingie[-1], tip_loc_y] 
        wl_zloc = [z_loc_wingie[-1],tip_loc_z] 
        wl_xloc = [x_loc_wingie[-1], tip_loc_x] 
#        print("hier")
#        print(wl_xloc, wl_yloc,wl_zloc)
         
        chord_wingie = tip_chord - (location_along_radius/winglet_b)*(tip_chord-winglet_tip) 
        wl_chord = [chord_wingie[-1], winglet_tip] 
#        print(wl_chord)
        return chord_wingie, wl_chord, wl_xloc, wl_yloc, wl_zloc
    
    chord_wingie, wl_chord, wl_xloc, wl_yloc, wl_zloc = winglet_tip(tip_chord,taper_winglet, A,x_loc_wingie,z_loc_wingie,y_loc_wingie,angle_w,winglet_LE_sweep, location_along_radius)
#    print(wl_xloc,wl_yloc,wl_zloc)
    
    
#    print(chord_wingie, wl_chord)
#    print(wl_xloc,wl_yloc,wl_zloc)
    def avl_wing(M_cruise, S, root_chord, span, CD_0, chordwise_discretize_point, spanwise_discretize_points, Angle, x_loc_LE, y_loc_LE, z_loc_LE,\
                 chords, Ainc,Nspanwise, Sspace, x_loc_wingie, y_loc_wingie, z_loc_wingie, chord_wingie, wl_xloc, wl_yloc, wl_zloc, wl_chord,x_engine, y_engine, z_engine, engine_chords):
        print("check")
        a = 3.34
        b = 8.26
        tail_length = 19.11
        fus_length = 49.3
        nosso_length = 8.3
        
        x = np.arange(0,3.34,0.5)
        
        y = -(np.sqrt((1-x**2/a**2)*b**2)-8.26)
        y_tail = tail_length-tail_length/a*x
        a = np.append(x, [3.34])
        b = np.append(y, [8.26])
        btail = np.append(y_tail, [0.0])  + (fus_length-tail_length)

        y_fus = a
        x_fus = b
        z_fus = np.full(len(y_fus),-3.34)
        chords_fus = btail-b
        
        with open("winglet.avl", "w") as text_file: 
            print("Wing" + "\n"
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
                  "Wing", "\n" +
                  str(chordwise_discretize_point), "1.0 " + str(spanwise_discretize_points), "-2.0" + "\n"
                  "YDUPLICATE" + "\n" +
                  str(0.0), "\n" +
                  "ANGLE" + "\n" +
                  str(Angle), file=text_file)
            for i in range(3):
                print("SECTION", file=text_file)
                print(round(x_loc_LE[i], 3), round(y_loc_LE[i], 3), round(z_loc_LE[i]-1, 3), round(chords[i], 3), Ainc[i],
                      Nspanwise[i], Sspace[i], file=text_file)
                print("CDCL" + "\n"
                      "-0.5 0.0069 0.3 0.0067 0.9 0.0069", file=text_file)
                print("CDCL" + "\n"
                      "-0.91 0.012 0.31 0.008 1.41 0.013", file=text_file)
                if i < 1:
                    print("AFILE" + "\n"
                                    "SC20616.dat", file=text_file)    
                else:
                    print("AFILE" + "\n"
                          "sc20612.dat", file=text_file)  
            
            print("SURFACE" + "\n" + 
                  "Winglet", "\n" + 
                  str(chordwise_discretize_point), "1.0 " + str(spanwise_discretize_points), "-2.0" + "\n" 
                  "COMPONENT", "\n" + 
                  str(1), "\n" 
                  "YDUPLICATE" + "\n" + 
                  str(0.0), "\n" + 
                  "ANGLE" + "\n" + 
                  str(Angle), file=text_file) 
            for i in range(np.size(x_loc_wingie)): 
                print("SECTION", file=text_file) 
                print(round(x_loc_wingie[i], 2), round(y_loc_wingie[i], 2), round(z_loc_wingie[i]-1, 2), round(chord_wingie[i],2), str(0.0),str(0.0),str(0.0), file=text_file) 
                   
            print("\n" + "SURFACE" + "\n" + 
              "Winglet2", "\n" + 
              str(chordwise_discretize_point), "1.0 " + str(spanwise_discretize_points), "-2.0" + "\n" 
              "COMPONENT", "\n" + 
              str(1), "\n" 
                  "YDUPLICATE" + "\n" + 
                  str(0.0), "\n" + 
              "ANGLE" + "\n" + 
              str(Angle), file=text_file) 
            for i in range(np.size(wl_xloc)): 
                print("SECTION", file=text_file) 
                print(round(wl_xloc[i], 2), round(wl_yloc[i], 2), round(wl_zloc[i]-1, 2), round(wl_chord[i],2), str(0.0),str(0.0),str(0.0), file=text_file) 
         
                print("CDCL" + "\n"
                      "-0.91 0.012 0.31 0.008 1.41 0.013", file=text_file) 
                print("AFILE" + "\n"
                      "SC(2)-0010.txt", file=text_file) 
#            print(wl_yloc)
##            
            print("SURFACE" + "\n" +
              "Stab", "\n" +
              str(chordwise_discretize_point), "1.0 " + str(spanwise_discretize_points), "-1.1" + "\n"
                                                                                                  "YDUPLICATE" + "\n" +
              str(0.0), "\n" +
              "ANGLE" + "\n" +
              str(0.0), file=text_file)
            for i in range(2):
                print("SECTION", file=text_file)
                print(round(x_loc_LE_h[i], 3), round(y_loc_LE_h[i], 3), round(z_loc_LE_h[i]-3.34, 3), round(h_chords[i], 3),
                      str(0.0), file=text_file)
            print("AFILE" + "\n"
                            "n0010.dat", file=text_file)

#            print("\n" + "SURFACE" + "\n" +
#                  "FIN", "\n" +
#                  str(chordwise_discretize_point), "1.0 " + str(spanwise_discretize_points), "1.0" + "\n"
#                                                                                                     "YDUPLICATE" + "\n" +
#                  str(0.0), "\n" +
#                  "ANGLE" + "\n" +
#                  str(0.0), file=text_file)
#            for i in range(2):
#                print("SECTION", file=text_file)
#                print(round(x_loc_LE_v[i], 3), round(y_loc_LE_v[i], 3), round(z_loc_LE_v[i]-3.34, 3), round(v_chords[i], 3),
#                      str(0.0), file=text_file)
#            print("AFILE" + "\n"
#                            "n0010.dat", file=text_file)
            print("\n" + "SURFACE" + "\n" +
                  "NACELLE", "\n" +
                  str(6.0), "1.0 " + str(12.0), "0.0" + "\n"
                  "YDUPLICATE" + "\n" +
                  str(0.0), "\n" +
                  "ANGLE" + "\n" +
                  str(0.0), file=text_file)   
            for i in range(len(x_engine)):
                print("SECTION", file=text_file)
                print(round(x_engine[i], 3), round(y_engine[i], 3), round(z_engine[i]-1, 3), round(engine_chords[i], 3),
                      str(0.0), str(1), str(0.), file=text_file)
            print("\n" + "SURFACE" + "\n" +
                  "Fuselage H" + "\n" +
                  str(24), str(1.0), "\n"+
                  "COMPONENT" +"\n" +
                  str(6), "\n" +
                  "YDUPLICATE" + "\n" +
                  str(0.0), file=text_file)
            for i in range(len(chords_fus)):
                 print("SECTION", file=text_file)
                 print(round(x_fus[i], 3), round(y_fus[i], 3), round(z_fus[i], 3), round(chords_fus[i], 3),
                      str(0.0), str(1), str(0.), file=text_file)
            
            print("\n" + "SURFACE" + "\n" +
                  "Fuselage top" + "\n" +
                  str(24), str(1.0), "\n"+
                  "NOWAKE" + "\n" + 
                  "COMPONENT" +"\n" +
                  str(6), file=text_file)
            for i in range(len(chords_fus)):
                 print("SECTION", file=text_file)
                 print(round(x_fus[i], 3), str(0.0), round(y_fus[i]-3.34, 3), round(chords_fus[i], 3),
                      str(0.0), str(1), str(0.), file=text_file)
            
            print("\n" + "SURFACE" + "\n" +
                  "Fuselage bottom" + "\n" +
                  str(24), str(1.0), "\n"+
                  "COMPONENT" +"\n" +
                  str(6), file=text_file)
            for i in range(len(chords_fus)):
                 print("SECTION", file=text_file)
                 print(round(x_fus[i], 3), str(0.0), round(-y_fus[i]-3.34, 3), round(chords_fus[i], 3),
                      str(0.0), str(1), str(0.), file=text_file)
            print("\n" + "SURFACE" + "\n" +
                  "STRUT"+ "\n" + 
                  str(24), str(1.0), "\n" +
                  "YDUPLICATE" + "\n" +
                  str(0.0), "\n" +
                  "COMPONENT" + "\n" +
                  str(5), file=text_file)
            
            for i in range(2):
                 print("SECTION", file=text_file)
                 print(round(x_strut[i], 3), round(y_strut[i],3), round(z_strut[i], 3), round(strut_chords[i], 3),
                      str(-1.86), str(1), str(0.), file=text_file)
                 print("AFILE" + "\n"
                "n0010.dat", file=text_file)
                 
    avl_wing(M_cruise, S, root_chord, span, CD_0, chordwise_discretize_point, spanwise_discretize_points, Angle, x_loc_LE, y_loc_LE, z_loc_LE,\
                 chords, Ainc,Nspanwise, Sspace, x_loc_wingie, y_loc_wingie, z_loc_wingie, chord_wingie, wl_xloc, wl_yloc, wl_zloc, wl_chord, x_engine, y_engine, z_engine, engine_chords) 
    
#avl(root_chord, tip_chord,span, LE_sweep, dihedral, 10,radius,R_points, winglet_LE_sweep, taper_winglet, A, intermed_chord)

def find_clalpha(M_cruise, CD_0, filename):
    define_run_condition(M_cruise, CD_0)
    ROOT_dir = os.path.dirname(os.path.abspath("winglet.py"))
#    print(ROOT_dir)
    alpha_range = [0, 10]
    CL_range = []
    for j in range(len(alpha_range)):
        p = subprocess.Popen(str(ROOT_dir) + "/avl.exe", stdin=subprocess.PIPE, stdout=subprocess.DEVNULL,
                             universal_newlines=True)
        set_alpha = "a a " + str(alpha_range[j])
        p.communicate(os.linesep.join(["load", filename, "case", "mach" + str(M_cruise) + ".run", "oper", set_alpha, "x", "ft", "endresult"]))
        lines = [line.rstrip('\n') for line in open('endresult')]
        CL = float(lines[23].split()[2])
        CL_range.append(CL)
        os.remove("endresult")
    CL_alpha = round((CL_range[1] - CL_range[0]) / (alpha_range[1] - alpha_range[0]), 5)
#    os.remove("mach" + str(M_cruise) + ".run")
    return CL_alpha

a = find_clalpha(0.0,CD_0,"winglet")
print("alpha", a, a*57.3)

def run_avl(cl_cruise, M_cruise, CD_0):
#    define_run_condition(M_cruise, CD_0)    
    ROOT_dir = os.path.dirname(os.path.abspath("winglet.py"))
    p = subprocess.Popen(str(ROOT_dir) + "/avl.exe", stdin=subprocess.PIPE, stdout=subprocess.DEVNULL,
                        universal_newlines=True)
    set_cl_cruise = "a c " + str(cl_cruise)
    p.communicate(os.linesep.join(
        ["load", "winglet", "case", "mach0.75.run", "oper", set_cl_cruise, "x", "ft", "endresult"]))
    lines = [line.rstrip('\n') for line in open('endresult')]
    alpha = float(lines[15].split()[2])
    CD = float(lines[24].split()[2])
    CDi = float(lines[25].split()[5])
    print(CDi)
    print(CD)
    print(alpha)
#    e = float(lines[27].split()[5])
    # print("Alpha is ", alpha)
    # print("CD is ", CD)
    # print("eff factor is ", e)
    os.remove("endresult")
#    os.remove("mach" + str(M_cruise) + ".run")
    return CD

#
angle_w = np.arange(10,20,10)
cla = []
cdo = []
#step = 0.2
c = np.arange(.62,.72,.1)
for i in angle_w:
    print(i)
    engine_data = engine_avl(engine_placement)
    avl(root_chord, tip_chord,span, LE_sweep, dihedral, i,radius,R_points, winglet_LE_sweep, taper_winglet, A, intermed_chord, engine_data[0], engine_data[1], engine_data[2], engine_data[3])
    cla.append( find_clalpha(0, CD_0, "winglet.avl"))
    for j in c:
        print(j)
        cdo.append(run_avl(j,0,CD_0))
##%%
#x = np.reshape(np.asarray(cdo), (-1,5))        
#
#for reeks in x:
#    plt.figure(2)
#    plt.plot(reeks,c)
#plt.legend(['0 deg','10 deg','20 deg','30 deg','40 deg'])
##plt.ylim((1.5,2.5))
##plt.xlim((0.006,0.012))
#plt.xlabel("$CD_i$ [-]")
#plt.ylabel(r'$\alpha$ [deg]')
#        plt.plot(a,b)
#        plt.plot((3.34,3.34), (nosso_length,(fus_length-tail_length)))
#        plt.plot(a, btail)
#        plt.axis('equal')


#
#from matplotlib import pyplot as plt
#step = .1
#CL = np.arange(0,1.5+step,step)
#CD = np.zeros(len(CL))
#for i in range(len(CL)):
#    print(i)
##    define_run_condition(M_cruise,CD_0)
#    p = subprocess.Popen(r'C:\Users\floyd\OneDrive\Documenten\GitHub\Bigger_is_Better\avl\avl.exe', stdin=subprocess.PIPE, universal_newlines=True)
#    set_CL = "a c " + str(CL[i])
#    p.communicate(os.linesep.join(["load", "winglet","case","mach0.75.run", "oper", set_CL, "x","ft", "endresult"]))          
#    
#    lines = [line.rstrip('\n') for line in open('endresult')]
#    CD[i] = float(lines[24].split()[2])
#    CDi = float(lines[25].split()[5])
#    print(CDi)
#    os.remove("endresult")
#    print("check")
##    os.remove("mach" + str(M_cruise) + ".run")
#plt.figure(2)
#plt.plot(CD,CL)
#CD2 = 0.0233 + CL**2/(np.pi*13.82*0.72)
##plt.plot(CD2,CL)
#plt.grid()    
#plt.xlabel("CD [-]")
#plt.ylabel("CL [-]")
##plt.legend(["AVL","Theoretical"])
#plt.ylim((0,1.6))
#plt.show()        