import numpy as np
import constants_and_conversions as cc
import matplotlib.pyplot as plt

x_strut = 20
x_eng = 7
d_fus = 5.95
b = 30
tc = 0.129

E_wing = 69 * (10 ** 9)
I_wing = 0.0025
	
gamma = np.arctan(d_fus / x_strut)
#print(np.degrees(gamma), np.sin(gamma))

A_strut = 0.25 * np.pi * ((4.5 / 100) ** 2)
L_strut = np.sqrt(d_fus ** 2 + x_strut ** 2)
E_strut = 69 * (10 ** 9)

l_max = 42580
l_min = 8011
l_slope = (l_min - l_max) / b

w_max = 986
w_min = 95
w_slope = (w_min - w_max) / b

w_eng = 73144
w_fuel = 200716 / 2
w_wing = 205725

length = b
thickness = 5 / 1000
Cr = 6.14
Ct = 1.83

tot_volume = (0.5 * ((Cr ** 2) * tc + (Ct ** 2) * tc) * (b / 2)) * 2  # m^3
w_spec_w = w_wing / tot_volume


def calc_chord(cr, ct, b, x):
    return cr - ((cr - ct) / (b / 2)) * x


def calc_volume(tc, chord, width):
    return chord * tc * chord * width * 0.5


def calc_fuel_range(start, tc, dx):
    V_fuel = (w_fuel / cc.g_0) / (0.804 * 1000)
    forces = []
    front = []
    back = []
    
    while V_fuel > 0:        
        chord = calc_chord(Cr, Ct, b * 2, start - dx / 2)
        volume = calc_volume(tc, chord, dx)

        V_fuel -= volume
        
        forces.append(volume * (0.804 * 1000) * cc.g_0)
        if start < x_strut:
             front.append(volume * (0.804 * 1000) * cc.g_0)
        else:
            back.append(volume * (0.804 * 1000) * cc.g_0)
        start -= dx

#    print("End of fuel ", start)
#    print(np.sum(np.array(front)))
#    print(np.sum(np.array(back)))
    return start, forces[0], forces[-1]   


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


def find_strut_force():
    diff = 0.1 / 1000
    strut_force = 0
    
    chord = calc_chord(Cr, Ct, b * 2, x_strut)
    I_wing = MOI2(chord, thickness)
    start_fuel = 28
    
    chord = calc_chord(Cr, Ct, b * 2, start_fuel)
    volume = calc_volume(tc, chord, 1)
    
    f_min = volume * (0.804 * 1000) * cc.g_0
    V_remaining = ((w_fuel - f_min * start_fuel) / cc.g_0) / (0.804 * 1000)

    f_diff = V_remaining * cc.g_0 * (0.804 * 1000)
    f_max = f_min + (f_diff * 2) / start_fuel
    f_slope = (f_min - f_max) / start_fuel


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
        
        fuel_aft_cst = -1 * f_min * (start_fuel - x_strut)
        d_fuel_aft = fuel_aft_cst * (x_strut ** 3) / (3 * E_wing * I_wing)
        d_wing += d_fuel_aft
        
        mom_aft_cst = lift_aft_cst * (b - x_strut) / 2 + weight_aft_cst * (b - x_strut) / 2 + fuel_aft_cst * (start_fuel - x_strut) / 2
        d_mom_aft = mom_aft_cst * (x_strut ** 2) / (2 * E_wing * I_wing)
        d_wing += d_mom_aft 
        
        a = x_strut - x_eng
        d_eng = -1 * w_eng * (2 * (x_strut ** 3) - 3 * (x_strut ** 2) * a + (a ** 3)) / (6 * E_wing * I_wing)
        d_wing += d_eng
        
        lift_front_cst = l_max + l_slope * x_strut 
        d_lift_front = lift_front_cst * (x_strut ** 4) / (8 * E_wing * I_wing)
        d_wing += d_lift_front 
        
        weight_front_cst = -1 * (w_max + w_slope * x_strut)
        d_weight_front = weight_front_cst * (x_strut ** 4) / (8 * E_wing * I_wing)
        d_wing += d_weight_front
        
        fuel_front_cst = -1 * (f_max + f_slope * x_strut) # (f_max + f_slope * (x_strut - end_fuel)) 
        d_fuel_front = fuel_front_cst * (x_strut ** 4) / (8 * E_wing * I_wing)
        d_wing += d_fuel_front
        
        lift_front_diff = l_max - lift_front_cst
        d_lift_front_diff = lift_front_diff * (x_strut ** 4) / (30 * E_wing * I_wing)
        d_wing += d_lift_front_diff
        
        weight_front_diff = -1 * (w_max + weight_front_cst)
        d_weight_front_diff = weight_front_diff * (x_strut ** 4) / (30 * E_wing * I_wing)
        d_wing += d_weight_front_diff
        
        fuel_front_diff = -1 * (f_max + fuel_front_cst)
        d_fuel_front_diff = fuel_front_diff * (x_strut ** 4) / (30 * E_wing * I_wing)
        d_wing += d_fuel_front_diff
        
        lift_aft_diff = 0.5 * (lift_front_cst - l_min) * (b - x_strut)
        d_lift_aft_diff = lift_aft_diff * (x_strut ** 3) / (3 * E_wing * I_wing)
        d_wing += d_lift_aft_diff
        
        weight_aft_diff = -0.5 * (-weight_front_cst - w_min) * (b - x_strut)
        d_weight_aft_diff = weight_aft_diff * (x_strut ** 3) / (3 * E_wing * I_wing)
        d_wing += d_weight_aft_diff
        
        fuel_aft_diff = -0.5 * (-fuel_front_cst - f_min) * (start_fuel - x_strut)
        d_fuel_aft_diff = fuel_aft_diff * (x_strut ** 3) / (3 * E_wing * I_wing)
        d_wing += d_fuel_aft_diff
        
        mom_aft_diff = lift_aft_diff * (b - x_strut) / 3 + weight_aft_diff * (b - x_strut) / 3 + fuel_aft_diff * (start_fuel - x_strut) / 3
        d_mom_aft_diff = mom_aft_diff * (x_strut ** 2) / (2 * E_wing * I_wing)
        d_wing += d_mom_aft_diff 
        
        d_strut = (strut_force * L_strut) / (E_strut * A_strut)
        
        if d_wing - d_strut < diff:
            break
    
    lift = lift_front_cst * x_strut + 0.5 * lift_front_diff * x_strut
    weight = weight_front_cst * x_strut + 0.5 * weight_front_diff * x_strut + weight_aft_cst + weight_aft_diff
    fuel_weight = fuel_front_cst * x_strut + 0.5 * fuel_front_diff * x_strut
    shear_strut = lift_aft_cst + weight_aft_cst + lift_aft_diff + weight_aft_diff + fuel_aft_cst + fuel_aft_diff
    
    V_root_y = lift + weight + fuel_weight - w_eng + shear_strut - strut_force * np.sin(gamma) 
        
    mom_aft_strut = -1 * (mom_aft_cst + mom_aft_diff)
    mom_shear_aft_strut = -1 * shear_strut * x_strut
    mom_eng = x_eng * w_eng
    mom_lift = -lift_front_cst * (x_strut ** 2) / 2 - 0.5 * lift_front_diff * (x_strut ** 2) / 3
    mom_weight = -1 * weight_front_cst * (x_strut ** 2) / 2 + 0.5 * weight_front_diff * (x_strut ** 2) / 3
    mom_fuel = -1 * fuel_front_cst * (x_strut ** 2) / 2 + 0.5 * fuel_front_diff * (x_strut ** 2) / 3
    mom_strut_force = np.sin(gamma) * strut_force * x_strut
    
    M_root = -1 * (mom_aft_strut + mom_shear_aft_strut + mom_eng + mom_lift + mom_weight + mom_strut_force + mom_fuel)
    
    print("Strut force in N: ", strut_force)
    print("Wing deflection in mm: ", d_wing * 1000)
    print("Strut deflection in mm: ", d_strut * 1000)
    print("Strut stress in mpa: ", (strut_force / A_strut) / (10 ** 6))
    
    print()
    print("The root moment in N: ", M_root)
    print("The shear force at the root: ", V_root_y)


    return strut_force, f_min, f_max, start_fuel


def defl_strutted(strut_force, f_min, f_max, start_fuel):
    lift_cst = l_min
    lift_diff = l_max - l_min
    
    weight_cst = w_min
    weight_diff = w_max - w_min
    
    fuel_cst = f_min
    fuel_diff = f_max - f_min

    x_range = np.arange(0, length + 1 / 1000, 1 / 1000)
    deflections = []
    deflect = [[], [], [], [], [], [], []]

    for x in x_range:
        
        chord = calc_chord(Cr, Ct, b * 2, x)
        I_wing = MOI2(chord, thickness)
                
        # Lift input
        lft_cst = (lift_cst * (x ** 2) / (24 * E_wing * I_wing)) * ((x ** 2) - 4 * length * x + 6 * (length ** 2))
        deflect[0].append(lft_cst)
        defl = lft_cst
        
        lft_diff = (lift_diff * (x ** 2) / (120 * E_wing * I_wing * length)) * (10 * (length ** 3) - 10 * (length ** 2) * x + 5 * length * (x ** 2) - x ** 3)
        deflect[1].append(lft_diff)
        defl += lft_diff
        
        # Weight input
        w_cst = ((-weight_cst * (x ** 2)) / (24 * E_wing * I_wing)) * ((x ** 2) - 4 * length * x + 6 * (length ** 2))
        deflect[2].append(w_cst)
        defl += w_cst
        
        w_diff = ((-weight_diff * (x ** 2)) / (120 * E_wing * I_wing * length)) * (10 * (length ** 3) - 10 * (length ** 2) * x + 5 * length * (x ** 2) - x ** 3)
        deflect[3].append(w_diff)
        defl += w_diff
        
        # Fuel input
        a_fuel = length - start_fuel
        y_fuel_cst = (-fuel_cst / (24 * E_wing * I_wing)) * ((length - a_fuel) ** 3) * (3 * length + a_fuel)
        theta_fuel_cst = (fuel_cst / (6 * E_wing * I_wing)) * ((length - a_fuel) ** 3) * (length - x)
        if length - x > a_fuel:
            mac_fuel_cst = ((-fuel_cst) / (24 * E_wing * I_wing)) * ((length - x - a_fuel ) ** 4)
        else:
            mac_fuel_cst = 0
        
        f_cst = y_fuel_cst + theta_fuel_cst + mac_fuel_cst
        
        y_fuel_diff = (-fuel_diff / (120 * E_wing * I_wing)) * ((length - a_fuel) ** 3) * (4 * length + a_fuel)
        theta_fuel_diff = (fuel_diff / (24 * E_wing * I_wing)) * ((length - a_fuel) ** 3) * (length - x)
        if length - x > a_fuel:
            mac_fuel_diff = (-fuel_diff / (120 * E_wing * I_wing * (length - a_fuel))) * ((length - x - a_fuel ) ** 5)
        else:
            mac_fuel_diff = 0
        
        f_diff = y_fuel_diff + theta_fuel_diff + mac_fuel_diff
        f_defl = f_cst + f_diff
        deflect[4].append(f_defl)
        defl += f_defl
        
        # Strut input
        a_str = length - x_strut
        y_a_str = -strut_force * np.sin(gamma) * (2 * (length ** 3) - 3 * (length ** 2) * a_str + (a_str ** 3)) / (6 * E_wing * I_wing)
        theta_a_str = ((strut_force * np.sin(gamma) * ((length - a_str) ** 2)) / (2 * E_wing * I_wing)) * (length - x)
        
        if length - x > a_str:
            mac_a_str = (-(strut_force * np.sin(gamma)) / (6 * E_wing * I_wing)) * ((length - x - a_str) ** 3)
        else:
            mac_a_str = 0
            
        diff_str = y_a_str + theta_a_str + mac_a_str
        deflect[5].append(diff_str)
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
        deflect[6].append(diff_eng)
        defl += diff_eng
        
        deflections.append(defl)
        
    return x_range, deflections, deflect


f_str, f_min, f_max, start = find_strut_force()
#x_range, defl_non_str, deflects_non = defl_strutted(0, f_min, f_max, start)
x_range, defl_str, deflects_str = defl_strutted(f_str, f_min, f_max, start)
labels = ["lift cst", "lift diff", "w cst", "w diff", "fuel bot", "fuel top", "strut", "engine"]

for i in range(len(deflects_str)):
    plt.plot(x_range, deflects_str[i], label= labels[i])
#    plt.plot(x_range, deflects_non[i], label= labels[i])
    
plt.plot(x_range, defl_str)
#plt.plot(x_range, defl_non_str)
#plt.axis("scaled")
plt.legend()
plt.xlabel("Spanwise location [m]")
plt.ylabel("Deflection [m]")
plt.show()

print()
print(defl_str[20000] * 1000)
#print(defl_non_str[20000] * 1000)


