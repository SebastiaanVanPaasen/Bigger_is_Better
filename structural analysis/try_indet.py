import numpy as np
import matplotlib.pyplot as plt

x_strut = 15
x_eng = 7
d_fus = 6.8
b = 25

E_wing = 69 * (10 ** 9)
I_wing = 0.005
	
gamma = np.arctan(d_fus / x_strut)
#print(np.degrees(gamma), np.sin(gamma))

A_strut = 0.25 * np.pi * ((4.1 / 100) ** 2)
L_strut = np.sqrt(d_fus ** 2 + x_strut ** 2)
E_strut = 69 * (10 ** 9)

l_max = 42580
l_min = 8011
l_slope = (l_min - l_max) / b

w_max = 986
w_min = 95
w_slope = (w_min - w_max) / b

w_eng = 73144

length = b


def defl_strutted():
    diff = 0.00001
    strut_force = 0
    
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
        
        mom_aft_cst = lift_aft_cst * (b - x_strut) / 2 + weight_aft_cst * (b - x_strut) / 2
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
        
        lift_front_diff = l_max - lift_front_cst
        d_lift_front_diff = lift_front_diff * (x_strut ** 4) / (30 * E_wing * I_wing)
        d_wing += d_lift_front_diff
        
        weight_front_diff = w_max + weight_front_cst
        d_weight_front_diff = -1 * weight_front_diff * (x_strut ** 4) / (30 * E_wing * I_wing)
        d_wing += d_weight_front_diff
        
        lift_aft_diff = 0.5 * (lift_front_cst - l_min) * (b - x_strut)
        d_lift_aft_diff = lift_aft_diff * (x_strut ** 3) / (3 * E_wing * I_wing)
        d_wing += d_lift_aft_diff
        
        weight_aft_diff = 0.5 * (weight_front_cst + w_min) * (b - x_strut)
        d_weight_aft_diff = weight_aft_diff * (x_strut ** 3) / (3 * E_wing * I_wing)
        d_wing += d_weight_aft_diff
        
    
        mom_aft_diff = lift_aft_diff * (b - x_strut) / 3 + weight_aft_diff * (b - x_strut) / 3
        d_mom_aft_diff = mom_aft_diff * (x_strut ** 2) / (2 * E_wing * I_wing)
        d_wing += d_mom_aft_diff  
        
        
        d_strut = (strut_force * L_strut) / (E_strut * A_strut)
        
        if d_wing - d_strut < diff:
            break
    
    
    lift = lift_front_cst * x_strut + 0.5 * lift_front_diff * x_strut
    weight = weight_front_cst * x_strut + 0.5 * weight_front_diff * x_strut
    shear_strut = lift_aft_cst + weight_aft_cst + lift_aft_diff + weight_aft_diff
    
    V_root_y = lift + weight + shear_strut - strut_force * np.sin(gamma)
    
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
    print("Strut stress in mpa: ", (strut_force / A_strut) / (10 ** 6))
    
    x_range = np.arange(0, length, 1 / 1000)
    deflections = []
    for x in x_range:
        
        defl = (-strut_force * np.sin(gamma) * (x ** 2) / (6 * E_wing * I_wing)) * (3 * length - x)
        defl += -mom_aft_strut * (x ** 2) / (2 * E_wing * I_wing)
        defl += (shear_strut * (x ** 2) / (6 * E_wing * I_wing)) * (3 * length - x)
        defl += (lift_front_cst * (x ** 2) / (24 * E_wing * I_wing)) * ((x ** 2) - 4 * length * x + 6 * (length ** 2))
        defl += (lift_front_diff * (x ** 2) / (120 * E_wing * I_wing * length)) * (10 * (length ** 3) - 10 * (length ** 2) * x + 5 * length * (x ** 2) - x ** 3)
        defl += (weight_front_cst * (x ** 2) / (24 * E_wing * I_wing)) * ((x ** 2) - 4 * length * x + 6 * (length ** 2))
        defl += (-weight_front_diff * (x ** 2) / (120 * E_wing * I_wing * length)) * (10 * (length ** 3) - 10 * (length ** 2) * x + 5 * length * (x ** 2) - x ** 3)
    
        y_a = d_eng
        theta_a = ((w_eng * ((length - a) ** 2)) / (2 * E_wing * I_wing)) * (length - x)
        
        if length - x > a:
            mac_a = (-w_eng / (6 * E_wing * I_wing)) * ((length - x - a) ** 3)
        else:
            mac_a = 0
            
        diff = y_a + theta_a + mac_a
        
        defl += diff
        
        deflections.append(defl * 1000)
        
    return x_range, deflections


def defl_non_strutted():
    lift_front_cst = l_max + l_slope * x_strut 
    lift_front_diff = l_max - lift_front_cst
    
    lift_aft_cst = l_min * (b - x_strut)
    lift_aft_diff = 0.5 * (lift_front_cst - l_min) * (b - x_strut)
    
    weight_front_cst = -1 * (w_max + w_slope * x_strut + w_min)
    weight_front_diff = w_max + weight_front_cst
    
    weight_aft_cst = -1 * w_min * (b - x_strut)
    weight_aft_diff = 0.5 * (weight_front_cst + w_min) * (b - x_strut)
    
    mom_aft_diff = lift_aft_diff * (b - x_strut) / 3 + weight_aft_diff * (b - x_strut) / 3
    mom_aft_cst = lift_aft_cst * (b - x_strut) / 2 + weight_aft_cst * (b - x_strut) / 2
    mom_aft_strut = -1 * (mom_aft_cst + mom_aft_diff)
    
    a = x_strut - x_eng
    d_eng = -1 * w_eng * (2 * (x_strut ** 3) - 3 * (x_strut ** 2) * a + (a ** 3)) / (6 * E_wing * I_wing)
    
    shear_strut = lift_aft_cst + weight_aft_cst + lift_aft_diff + weight_aft_diff
    
    x_range = np.arange(0, length, 1/ 1000)
    deflections = []
    for x in x_range:
        
        defl = -mom_aft_strut * (x ** 2) / (2 * E_wing * I_wing)
        defl += (shear_strut * (x ** 2) / (6 * E_wing * I_wing)) * (3 * length - x)
        defl += (lift_front_cst * (x ** 2) / (24 * E_wing * I_wing)) * ((x ** 2) - 4 * length * x + 6 * (length ** 2))
        defl += (lift_front_diff * (x ** 2) / (120 * E_wing * I_wing * length)) * (10 * (length ** 3) - 10 * (length ** 2) * x + 5 * length * (x ** 2) - x ** 3)
        defl += (weight_front_cst * (x ** 2) / (24 * E_wing * I_wing)) * ((x ** 2) - 4 * length * x + 6 * (length ** 2))
        defl += (-weight_front_diff * (x ** 2) / (120 * E_wing * I_wing * length)) * (10 * (length ** 3) - 10 * (length ** 2) * x + 5 * length * (x ** 2) - x ** 3)
    
        y_a = d_eng
        theta_a = ((w_eng * ((length - a) ** 2)) / (2 * E_wing * I_wing)) * (length - x)
        
        if length - x > a:
            mac_a = (-w_eng / (6 * E_wing * I_wing)) * ((length - x - a) ** 3)
        else:
            mac_a = 0
            
        diff = y_a + theta_a + mac_a
        
        defl += diff
        
        deflections.append(defl * 1000)
        
    return x_range, deflections


x_range, defl_str = defl_strutted()
x_range, defl_non_str = defl_non_strutted()
    
plt.plot(x_range, defl_str)
plt.plot(x_range, defl_non_str)
plt.xlabel("Spanwise location [m]")
plt.ylabel("Deflection [mm]")
plt.show()
#print(defl * 1000)

