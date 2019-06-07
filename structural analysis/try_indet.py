import numpy as np
import constants_and_conversions as cc
import matplotlib.pyplot as plt

x_strut = 20
x_eng = 7
d_fus = 6.8
b = 30
tc = 0.14

E_wing = 69 * (10 ** 9)
I_wing = 0.0025
	
gamma = np.arctan(d_fus / x_strut)
#print(np.degrees(gamma), np.sin(gamma))

A_strut = 0.25 * np.pi * ((5.5 / 100) ** 2)
L_strut = np.sqrt(d_fus ** 2 + x_strut ** 2)
E_strut = 69 * (10 ** 9)

l_max = 42580
l_min = 8011
l_slope = (l_min - l_max) / b

w_max = 986
w_min = 95
w_slope = (w_min - w_max) / b

w_eng = 73144
w_fuel = 207500 / 2
#print("tot fuel", w_fuel)

length = b
thickness = 5 / 1000
Cr = 6.37
Ct = 1.89


def calc_chord(cr, ct, b, x):
    return cr - ((cr - ct) / (b / 2)) * x


def calc_volume(tc, chord, width):
    return chord * tc * chord * width * 0.5


def calc_fuel_range(start, tc, dx):
    V_fuel = (w_fuel / cc.g_0) / (0.804 * 1000)
    forces = []
    
    while V_fuel > 0:
        chord = calc_chord(Cr, Ct, b * 2, start - dx / 2)
        volume = calc_volume(tc, chord, dx)

        V_fuel -= volume
        
        forces.append(volume * (0.804 * 1000) * cc.g_0)
        start -= dx

#    print("End of fuel ", start)
    return start, forces      
    
    
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


def find_strut_force():
    diff = 0.1 / 1000
    strut_force = 0
    
    chord = calc_chord(Cr, Ct, b * 2, x_strut)
    I_wing = MOI2(chord, thickness)
    start_fuel = 28
    end_fuel, fuel_forces = calc_fuel_range(start_fuel, tc, 1/1000)
    f_slope = (fuel_forces[0] - fuel_forces[-1]) / (start_fuel - end_fuel)
    
    while True:
        strut_force += 10
        d_wing = 0
        
        d_strut_v = (-1 * np.sin(gamma) * strut_force * (x_strut ** 3)) / (3 * E_wing * I_wing)
        d_wing += d_strut_v
        
        lift_aft_cst = l_min * (b - x_strut)
        d_lift_aft = lift_aft_cst * (x_strut ** 3) / (3 * E_wing * I_wing)
        d_wing += d_lift_aft
        
        weight_aft_cst = -1 * w_min * (b - x_strut)
        d_weight_aft = weight_aft_cst * (x_strut ** 3) / (3 * E_wing * I_wing)
        d_wing += d_weight_aft
        
#        fuel_aft_cst = -1 * fuel_forces[0] *  ((start_fuel - x_strut) * 1000)
#        d_fuel_aft = fuel_aft_cst * (x_strut ** 3) / (3 * E_wing * I_wing)
#        d_wing += d_fuel_aft
        
        mom_aft_cst = lift_aft_cst * (b - x_strut) / 2 + weight_aft_cst * (b - x_strut) / 2# + fuel_aft_cst * (start_fuel - x_strut) / 2
        d_mom_aft = mom_aft_cst * (x_strut ** 2) / (2 * E_wing * I_wing)
        d_wing += d_mom_aft 
        
        a = x_strut - x_eng
        d_eng = -1 * w_eng * (2 * (x_strut ** 3) - 3 * (x_strut ** 2) * a + (a ** 3)) / (6 * E_wing * I_wing)
        d_wing += d_eng
        
        lift_front_cst = l_max + l_slope * x_strut 
        d_lift_front = lift_front_cst * (x_strut ** 4) / (8 * E_wing * I_wing)
        d_wing += d_lift_front 
        
        weight_front_cst = -1 * (w_max + w_slope * x_strut + w_min)
        d_weight_front = weight_front_cst * (x_strut ** 4) / (8 * E_wing * I_wing)
        d_wing += d_weight_front
        
#        fuel_front_cst = -1 * (fuel_forces[-1] + f_slope * (x_strut - end_fuel)) * 1000 
#        d_fuel_front = fuel_front_cst * (x_strut ** 4) / (8 * E_wing * I_wing)
#        d_wing += d_fuel_front
        
        lift_front_diff = l_max - lift_front_cst
        d_lift_front_diff = lift_front_diff * (x_strut ** 4) / (30 * E_wing * I_wing)
        d_wing += d_lift_front_diff
        
        weight_front_diff = w_max + weight_front_cst
        d_weight_front_diff = -1 * weight_front_diff * (x_strut ** 4) / (30 * E_wing * I_wing)
        d_wing += d_weight_front_diff
        
#        fuel_front_diff = -1 * (fuel_forces[-1] * 1000 + fuel_front_cst)
#        d_fuel_front_diff = -1 * fuel_front_diff * (x_strut ** 4) / (30 * E_wing * I_wing)
#        d_wing += d_fuel_front_diff
        
        lift_aft_diff = 0.5 * (lift_front_cst - l_min) * (b - x_strut)
        d_lift_aft_diff = lift_aft_diff * (x_strut ** 3) / (3 * E_wing * I_wing)
        d_wing += d_lift_aft_diff
        
        weight_aft_diff = 0.5 * (weight_front_cst + w_min) * (b - x_strut)
        d_weight_aft_diff = weight_aft_diff * (x_strut ** 3) / (3 * E_wing * I_wing)
        d_wing += d_weight_aft_diff
        
#        fuel_aft_diff = 0.5 * (fuel_front_cst + fuel_forces[0] * 1000) * (start_fuel - x_strut)
#        d_fuel_aft_diff = fuel_aft_diff * (x_strut ** 3) / (3 * E_wing * I_wing)
#        d_wing += d_fuel_aft_diff
        
        mom_aft_diff = lift_aft_diff * (b - x_strut) / 3 + weight_aft_diff * (b - x_strut) / 3 #+ fuel_aft_diff * (start_fuel - x_strut) / 3
        d_mom_aft_diff = mom_aft_diff * (x_strut ** 2) / (2 * E_wing * I_wing)
        d_wing += d_mom_aft_diff  
        
        d_strut = (strut_force * L_strut) / (E_strut * A_strut)
        
        if d_wing - d_strut < diff:
            break
    
    lift = lift_front_cst * x_strut + 0.5 * lift_front_diff * x_strut
    weight = weight_front_cst * x_strut + 0.5 * weight_front_diff * x_strut
    shear_strut = lift_aft_cst + weight_aft_cst + lift_aft_diff + weight_aft_diff #+ fuel_aft_cst + fuel_aft_diff
#    fuel_weight = fuel_front_cst * (x_strut - end_fuel) + 0.5 * fuel_front_diff * (x_strut - end_fuel)
#    tot_fuel_weight = fuel_aft_cst + fuel_aft_diff + fuel_weight

    V_root_y = lift + weight + shear_strut - strut_force * np.sin(gamma) - w_eng
        
    mom_aft_strut = -1 * (mom_aft_cst + mom_aft_diff)
    mom_shear_aft_strut = -1 * shear_strut * x_strut
    mom_eng = x_eng * w_eng
    mom_lift = -lift_front_cst * (x_strut ** 2) / 2 - 0.5 * lift_front_diff * (x_strut ** 2) / 3
    mom_weight = -1 * weight_front_cst * (x_strut ** 2) / 2 + 0.5 * weight_front_diff * (x_strut ** 2) / 3
    mom_strut_force = np.sin(gamma) * strut_force * x_strut
    
    M_root = -1 * (mom_aft_strut + mom_shear_aft_strut + mom_eng + mom_lift + mom_weight + mom_strut_force)
    
    print("Strut force in N: ", strut_force)
    print("Wing deflection in mm: ", d_wing * 1000)
    print("Strut deflection in mm: ", d_strut * 1000)
    print("The root moment in N: ", M_root)
    print("The shear force at the root: ", V_root_y)
    print("Strut stress in mpa: ", (strut_force / A_strut) / (10 ** 6))
    
    return strut_force, fuel_forces, start_fuel, end_fuel, f_slope


def defl_strutted(strut_force, fuel, start_fuel, end_fuel, f_slope):
    
    lift_cst = l_min
    lift_diff = l_max - l_min
    weight_cst = w_min
    weight_diff = w_max - w_min
    fuel_cst = fuel[0] * 1000
    fuel_diff = (fuel[-1] - fuel[0]) * 1000
    fuel_cst_top = (fuel[-1] + f_slope * (x_strut - end_fuel)) * 1000

    x_range = np.arange(0, length, 1 / 1000)
    deflections = []
    for x in x_range:
        
        chord = calc_chord(Cr, Ct, b * 2, x)
        I_wing = MOI2(chord, thickness)
                
        defl = (lift_cst * (x ** 2) / (24 * E_wing * I_wing)) * ((x ** 2) - 4 * length * x + 6 * (length ** 2))
        defl += (lift_diff * (x ** 2) / (120 * E_wing * I_wing * length)) * (10 * (length ** 3) - 10 * (length ** 2) * x + 5 * length * (x ** 2) - x ** 3)
        defl += (-weight_cst * (x ** 2) / (24 * E_wing * I_wing)) * ((x ** 2) - 4 * length * x + 6 * (length ** 2))
        defl += (-weight_diff * (x ** 2) / (120 * E_wing * I_wing * length)) * (10 * (length ** 3) - 10 * (length ** 2) * x + 5 * length * (x ** 2) - x ** 3)

        # Fuel input
        a_fuel_bottom = length - start_fuel
        
        # Constant load acting downwards, equal to fuel minimum
        y_a_fuel_cst_bot = (-fuel_cst / (24 * E_wing * I_wing)) *((length - a_fuel_bottom) ** 3) * (3 * length + a_fuel_bottom) 
        theta_a_fuel_cst_bot = (fuel_cst / (6 * E_wing * I_wing)) * ((length - a_fuel_bottom) ** 3) * (length - x)
        
        if length - x > a_fuel_bottom:
            mac_a_fuel_cst_bot = (-fuel_cst / (24 * E_wing * I_wing)) * ((x - a_fuel_bottom) ** 4)
        else:
            mac_a_fuel_cst_bot = 0
            
        diff_fuel_cst_bot = y_a_fuel_cst_bot + theta_a_fuel_cst_bot + mac_a_fuel_cst_bot
        
        # Varying load acting downwards, equal to difference between fuel min and fuel max
        y_a_fuel_diff_bot = (-fuel_cst / (24 * E_wing * I_wing)) *((length - a_fuel_bottom) ** 3) * (3 * length + a_fuel_bottom) - (fuel_diff / (120 * E_wing * I_wing)) * ((length - a_fuel_bottom) ** 3) * (4 * length + a_fuel_bottom)
        theta_a_fuel_diff_bot = (fuel_cst / (6 * E_wing * I_wing)) * ((length - a_fuel_bottom) ** 3) * (length - x) + (fuel_diff / (24 * E_wing * I_wing)) * ((length - a_fuel_bottom) ** 3)
        
        if length - x > a_fuel_bottom:
            mac_a_fuel_diff_bot = (-fuel_cst / (24 * E_wing * I_wing)) * ((x - a_fuel_bottom) ** 4) - (fuel_diff / (120 * E_wing * I_wing * (length - a_fuel_bottom))) * ((length - a_fuel_bottom) ** 5)
        else:
            mac_a_fuel_diff_bot = 0
            
        diff_fuel_diff_bot = y_a_fuel_diff_bot + theta_a_fuel_diff_bot + mac_a_fuel_diff_bot
        
        # Constant counter acting load on top
        y_a_fuel_cst_top = (fuel_cst_top / (24 * E_wing * I_wing)) *((length - a_fuel_bottom) ** 3) * (3 * length + a_fuel_bottom) 
        theta_a_fuel_cst_top = (fuel_cst_top / (6 * E_wing * I_wing)) * ((length - a_fuel_bottom) ** 3) * (length - x)
        
        if length - x > a_fuel_bottom:
            mac_a_fuel_cst_bot = (-fuel_cst / (24 * E_wing * I_wing)) * ((x - a_fuel_bottom) ** 4)
        else:
            mac_a_fuel_cst_bot = 0
            
        diff_fuel_cst_bot = y_a_fuel_cst_bot + theta_a_fuel_cst_bot + mac_a_fuel_cst_bot
        
        # Varying load acting downwards, equal to difference between fuel min and fuel max
        y_a_fuel_diff_bot = (-fuel_cst / (24 * E_wing * I_wing)) *((length - a_fuel_bottom) ** 3) * (3 * length + a_fuel_bottom) - (fuel_diff / (120 * E_wing * I_wing)) * ((length - a_fuel_bottom) ** 3) * (4 * length + a_fuel_bottom)
        theta_a_fuel_diff_bot = (fuel_cst / (6 * E_wing * I_wing)) * ((length - a_fuel_bottom) ** 3) * (length - x) + (fuel_diff / (24 * E_wing * I_wing)) * ((length - a_fuel_bottom) ** 3)
        
        if length - x > a_fuel_bottom:
            mac_a_fuel_diff_bot = (-fuel_cst / (24 * E_wing * I_wing)) * ((x - a_fuel_bottom) ** 4) - (fuel_diff / (120 * E_wing * I_wing * (length - a_fuel_bottom))) * ((length - a_fuel_bottom) ** 5)
        else:
            mac_a_fuel_diff_bot = 0
            
        diff_fuel_diff_bot = y_a_fuel_diff_bot + theta_a_fuel_diff_bot + mac_a_fuel_diff_bot


        
        # Strut input
        a_str = length - x_strut
        y_a_str = -strut_force * np.sin(gamma) * (2 * (length ** 3) - 3 * (length ** 2) * a_str + (a_str ** 3)) / (6 * E_wing * I_wing)
        theta_a_str = ((strut_force * np.sin(gamma) * ((length - a_str) ** 2)) / (2 * E_wing * I_wing)) * (length - x)
        
        if length - x > a_str:
            mac_a_str = (-(strut_force * np.sin(gamma)) / (6 * E_wing * I_wing)) * ((length - x - a_str) ** 3)
        else:
            mac_a_str = 0
            
        diff_str = y_a_str + theta_a_str + mac_a_str
#        print(diff_str)
        
        defl += diff_str
        
        # Engine input
        a_eng = length - x_eng
        y_a_eng = -1 * w_eng * (2 * (length ** 3) - 3 * (length ** 2) * a_eng + (a_eng ** 3)) / (6 * E_wing * I_wing)
        theta_a_eng = ((w_eng * ((length - a_eng) ** 2)) / (2 * E_wing * I_wing)) * (length - x)
        
        if length - x > a_eng:
            mac_a_eng = (-w_eng / (6 * E_wing * I_wing)) * ((length - x - a_eng) ** 3)
        else:
            mac_a_eng = 0
            
        diff_eng = y_a_eng + theta_a_eng + mac_a_eng
#        print(diff_eng)
#        
        defl += diff_eng
        
        deflections.append(defl)
        
    return x_range, deflections


f_str, fuel, start, end, slope = find_strut_force()
#x_range, defl_non_str = defl_strutted(0)
x_range, defl_str = defl_strutted(f_str, fuel, start, end, slope)
#
plt.plot(x_range, defl_str)
##plt.plot(x_range, defl_non_str)
##plt.axis("scaled")
plt.xlabel("Spanwise location [m]")
plt.ylabel("Deflection [m]")
plt.show()
##
#
#print(defl_str[22500])
#print(defl_non_str[22500])


