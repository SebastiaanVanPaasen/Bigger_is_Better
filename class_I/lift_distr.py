import numpy as np
import subprocess
import os


#from loading_and_moment_diagrams import c
#ROOT_DIR = os.path.dirname(os.path.abspath("structural analysis"))

#S = wing["S"]#300.#286.02#184.16#193.72#220.27
#span = wing["b"]#60.#47.83#39.56#41.76#55.53
##.#8.#8.5#9.#14.
#taper = wing["Taper"]#0.4#0.31#0.4#0.31
##Sweep0 = m.atan(m.tan(wing["Sweep"]) - 4 / AR * (-0.25 * (1 - taper) / (1 + taper))) # rad
#qc_sweep = wing["Sweep"]
#dihedral = 0
#Cr = wing["C_root"]#8.#8.54#7.11#6.63#6.06(2*S)/((1+taper)*span)
#MAC = Cr*(2/3)*((1+taper+taper**2)/(1+taper))
#Ct = Cr*taper

def make_avl_file():
    # B777 used as reference aircraft
    S, span, taper, qc_sweep = 238, 60, 0.297, 0.5133 
    Cr = 6.14
    Ct = Cr * taper
    
    chords = [Cr, Ct]
    MAC = Cr*(2/3)*((1+taper+taper**2)/(1+taper))
    dihedral = 0
    
    CD_0 = 0.0177 
    Angle = 0.0
    
    dx = 0.25*Cr + span/2*np.tan(qc_sweep) - 0.25*Ct
    dz = span/2*np.tan(dihedral)
    
    x_loc_LE = [0, dx]
    y_loc_LE = [0, span/2]
    z_loc_LE = [0, dz]
    
    Ainc = [0.0, 0.0]
    spanwise_discretize_points = 50   #If you go too high then your computer is dead
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
            print("AFILE" + "\n""n2414.dat.txt", file=text_file)
            
            
#make_avl_file()

def lift_distribution(CL):        
    p = subprocess.Popen(r"C:\Users\mathi\Documents\DSE\Bigger_is_Better\avl\avl.exe", stdin=subprocess.PIPE, universal_newlines=True)
#    p = subprocess.Popen(r"C:\Users\sebas\OneDrive\Documents\DSE\Bigger_is_Better\avl\avl.exe", stdin=subprocess.PIPE, universal_newlines=True)

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
    
#output_avl = lift_distribution(0.8)


def get_correct_data(output_avl):
    MAC = 4.373
    x_pos = []
    cl = []
    cl_2 = []
    cd = []
    for i in range(len(output_avl)):
        x_pos.append(output_avl[i][1])
        cl.append(output_avl[i][4])#/c(output_avl[i][1]))
        cl_2.append(output_avl[i][4]/MAC)
        cd.append(output_avl[i][8])
    x_pos = x_pos[len(x_pos):int(len(x_pos)/2)-1:-1] + x_pos[0:int(len(x_pos)/2)]
    cl_2 = cl_2[len(x_pos):int(len(x_pos)/2)-1:-1] + cl_2[0:int(len(x_pos)/2)]
    cd = cd[len(x_pos):int(len(x_pos)/2)-1:-1] + cd[0:int(len(x_pos)/2)]
#    plt.scatter(x_pos,cl)
#    plt.scatter(y_pos,cd)
#    plt.grid()
#    print(cl, cl_2)
    return x_pos,cl_2, cd 
    
#x_pos,cl,cd = get_correct_data(output_avl)
