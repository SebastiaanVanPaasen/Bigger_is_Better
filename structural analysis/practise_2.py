import numpy as np
import constants_and_conversions as cc
import matplotlib.pyplot as plt


d_fus = 5.95
b = 30 
tc = 0.129

E_wing = 69 * (10 ** 9)

l_max = 42580
l_min = 8011
l_slope = (l_min - l_max) / b

w_max = 986
w_min = 95
w_slope = (w_min - w_max) / b

w_eng = 73144
w_fuel = 200716 / 2
w_wing = 205725
	
Cr = 6.14
Ct = 1.83

tot_volume = (0.5 * ((Cr ** 2) * tc + (Ct ** 2) * tc) * (b / 2)) * 2  # m^3
w_spec_w = w_wing / tot_volume

length = b
thickness = 5 / 1000


def calc_chord(cr, ct, b, x):
    return cr - ((cr - ct) / (b / 2)) * x


def calc_volume(tc, chord, width):
    return chord * tc * chord * width * 0.5   


dx = 0.1
x = 0
volume = 0
while x < 1:
    chord = calc_chord(Cr, Ct, b * 2, x + dx / 2)
    volume += calc_volume(tc, chord, dx)
    
    x += dx
    
w_max = volume * w_spec_w

x = 29
volume = 0
while x < 30:
    chord = calc_chord(Cr, Ct, b * 2, x + dx / 2)
    volume += calc_volume(tc, chord, dx)
    
    x += dx
    
w_min = volume * w_spec_w

w_slope = (w_min - w_max) / b


def MOI2(chordlength, airfoil_t):
#    x_centroid = 0.42 * chordlength
    y_centroid = 1.54E-2 * chordlength
    total_airfoil_thickness = 0.14 * chordlength
    area = 0.1 * chordlength * airfoil_t
    y1 = total_airfoil_thickness - y_centroid    
    y2 = total_airfoil_thickness - y1
#    x1 = 0.5 * chordlength - x_centroid
    Ixx = (0.2 * chordlength * (airfoil_t ** 3)) / 12 + ((y1 ** 2) * area) +((y2 ** 2) * area)
#    Iyy = (0.1 * (chordlength ** 3) * airfoil_t) / 12 + 2*((x1 ** 2) * area)
#    Ixy = (x1 * y1 * area + x1 * y2 * area)
#    J = 0.1 * chordlength * (airfoil_t ** 3) * (16 / 3 - 3.36 * airfoil_t/ (0.1 * chordlength) * (1 - 1 / 12 * (airfoil_t ** 4) / ((0.1 * chordlength) ** 4)))
#    print("Ixx: ",Ixx,"Iyy: " ,Iyy, "Ixy: ", Ixy,"J: ",J)
    return Ixx


def find_strut_force(x_strut, x_eng):
    A_strut = 0.25 * np.pi * ((4.5 / 100) ** 2)
    L_strut = np.sqrt(d_fus ** 2 + x_strut ** 2)
    E_strut = 69 * (10 ** 9)    
    gamma = np.arctan(d_fus / x_strut)

    diff = 0.1 / 1000
    strut_force = 0
    
    chord = calc_chord(Cr, Ct, b * 2, x_strut)
    I_wing = MOI2(chord, thickness)


    while True:
        strut_force += 10
        d_wing = 0
        
        d_strut_v = (-1 * np.sin(gamma) * strut_force * (x_strut ** 3)) / (3 * E_wing * I_wing)
        d_wing += d_strut_v
        
#        lift_front_cst = 100 
#        d_lift_front = lift_front_cst * (x_strut ** 4) / (8 * E_wing * I_wing)
#        d_wing += d_lift_front 
        
        lift_front_cst = 50
        d_lift_front = lift_front_cst * (x_strut ** 4) / (8 * E_wing * I_wing)
        d_wing += d_lift_front 
        
#        lift_aft_cst = 50 * (length - x_strut)
#        d_lift_aft = lift_aft_cst * (x_strut ** 3) / (3 * E_wing * I_wing)
#        d_wing += d_lift_aft
#        
#        mom_aft_cst = lift_aft_cst * 15 / 2 #+ weight_aft_cst * (b - x_strut) / 2 + fuel_aft_cst * (start_fuel - x_strut) / 2
#        d_mom_aft = mom_aft_cst * (x_strut ** 2) / (2 * E_wing * I_wing)
#        d_wing += d_mom_aft

        
        d_strut = (strut_force * L_strut) / (E_strut * A_strut)
        
        if d_wing - d_strut < diff:
            break
        
    print("Strut force in N: ", strut_force * np.sin(gamma))
#    print("Lift force in N: ", lift_front_cst * x_strut)
    print("Wing deflection in mm: ", d_wing * 1000)
    print("Strut deflection in mm: ", d_strut * 1000)
    print("Strut stress in mpa: ", (strut_force / A_strut) / (10 ** 6))
    
    return strut_force, gamma #, f_min, f_max, start_fuel


def defl_strutted(strut_force, gamma, x_strut, x_eng, f_min, f_max, start_fuel):
    lift_cst = 100
#    strut_force = 800
    lift = 0
#    dx = 1 / 1000
    x_range = np.arange(0, length + 1000 / 1000, 1000 / 1000)
    deflections = []
    deflect = [[], [], [], [], [], [], []]
    lift = []
    
    
    I_wing = MOI2(calc_chord(Cr, Ct, b*2, 30), thickness)
    
    a_lift = 14.5
    y_a_lift = (50 / (24 * E_wing * I_wing)) * ((length - a_lift) ** 3) *  (3 * length + a_lift)
    theta_a_lift = (-50 / (6 * E_wing * I_wing)) * ((length - a_lift) ** 3) 

    a_str = 15
    y_a_str = -strut_force * np.sin(gamma) * (2 * (length ** 3) - 3 * (length ** 2) * a_str + (a_str ** 3)) / (6 * E_wing * I_wing)
    theta_a_str = ((strut_force * np.sin(gamma) * ((length - a_str) ** 2)) / (2 * E_wing * I_wing)) 

#    a_lift_2 = 15
#    y_a_lift_2 = (50 / (24 * E_wing * I_wing)) * ((length - a_lift_2) ** 3) *  (3 * length + a_lift_2)
#    theta_a_lift_2 = (-50 / (6 * E_wing * I_wing)) * ((length - a_lift_2) ** 3) 
#    print(y_a_lift_2, theta_a_lift_2)
#    a_lift_3 = 15
#    y_a_lift_3 = (-50 / (24 * E_wing * I_wing)) * ((length - a_lift_3) ** 3) *  (3 * length + a_lift_3)
#    theta_a_lift_3 = (50 / (6 * E_wing * I_wing)) * ((length - a_lift_3) ** 3) 

    for x in x_range:
        
        chord = calc_chord(Cr, Ct, b * 2, x)
        I_wing = MOI2(chord, thickness)
                
        # Lift input      
        if x > a_lift:
            mac_lift = (50 / (24 * E_wing * I_wing)) * ((x - a_lift) ** 4)
        else:
            mac_lift = 0
            
        lft = y_a_lift + theta_a_lift * x + mac_lift       
        defl = lft
        
#        if x > a_lift_2:
#            mac_lift_2 = (50 / (24 * E_wing * I_wing)) * ((x - a_lift_2) ** 4)
##            print(mac_lift_2, x)
#        else:
#            mac_lift_2 = 0
#            
#            
#        lft_2 = y_a_lift_2 + theta_a_lift_2 * x + mac_lift_2
#        defl = lft_2
        
        
        
#        if x > a_lift:
#            mac_lift_3 = (-50 / (24 * E_wing * I_wing)) * ((x - a_lift_3) ** 4)
#        else:
#            mac_lift_3 = 0
#            
#        lft_3 = y_a_lift_3 + theta_a_lift_3 * x + mac_lift_3      
#        defl += lft_3
        
        

        
#        lft_diff = (lift_diff * (x ** 2) / (120 * E_wing * I_wing * length)) * (10 * (length ** 3) - 10 * (length ** 2) * x + 5 * length * (x ** 2) - x ** 3)
#        deflect[0].append(lft_cst + lft_diff)
                
#        deflect[0].append(lft)
#        deflect[1].append(lft_2)
#        deflect[2].append(theta_a_lift)
#        deflect[3].append(mac_lift)
        
        
        
        # Strut input
        a_str = 15
#        print(a_str)
        
        if x > 15:
#            print(x - a_str)
            mac_a_str = (-(strut_force * np.sin(gamma)) / (6 * E_wing * I_wing)) * ((x - a_str) ** 3)
        else:
            mac_a_str = 0
            
        diff_str = y_a_str + theta_a_str * x + mac_a_str
        
        deflect[2].append(diff_str)
        defl += diff_str
        
        deflections.append(defl)
        
    return x_range, deflections, deflect, lift


def find_deflection(x_strut, x_eng):
    f_str, angle = find_strut_force(x_strut, x_eng)
    x_range, defl_str, deflects_str, lift = defl_strutted(f_str, angle, x_strut, x_eng, 1, 1, 1)

#    print("Deflection according to equation in mm: ", defl_str[x_strut * 1000] * 1000)
    d = []
    for i in range(len(defl_str)):
        d.append(defl_str[len(defl_str) - 1 - i])

    print("deflection at the end", d[-1])

    plt.plot(x_range, d, label="tot")
#    plt.plot(x_new, lift)
#    plt.plot(x_new, deflects_str[0], label = "1")
#    plt.plot(x_new, deflects_str[1], label="2")
#    plt.plot(x_new, deflects_str[2], label="strut")
#    plt.plot(x_new, deflects_str[3], label="mac")
    plt.legend()
    plt.show()
    
#    return defl_str, x_range, deflects_str

find_deflection(15, 7)
#strut_range = [10, 15, 20, 22, 24, 26]
#engine_range = [6, 8, 10, 12, 14, 16]
#labels = ["lift", "weight", "fuel", "strut", "engine"]
#
#k = 1
#for i in range(len(strut_range)): 
#    plt.figure(i + 1)
#    for j in range(len(engine_range)):
#        
#        print()
#        print("Strut location: ", strut_range[i], "Engine location: ", engine_range[j])
#        print()
#        
#        defl, x_axis, components = find_deflection(strut_range[i], engine_range[j])
#        
##        for k in range(len(components)):
##            plt.plot(x_axis, components[k], label= labels[k])
#            
#        #    plt.plot(x_range, deflects_non[i], label= labels[i])
#            
#        plt.plot(x_axis, defl, label = "x_eng = " + str(engine_range[j]))        
#        #plt.plot(x_range, defl_non_str)
##        plt.axis("scaled")
#        
#        plt.title("strut at " + str(strut_range[i]))
#        plt.legend()
#        plt.xlabel("Spanwise location [m]")
#        plt.ylabel("Deflection [m]")
#        plt.show()
#    plt.savefig("strut position_" + str(strut_range[i]))
        
    

#x_range, defl_non_str, deflects_non = defl_strutted(0, f_min, f_max, start)
#print(defl_non_str[20000] * 1000)

