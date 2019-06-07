import numpy as np
import matplotlib.pyplot as plt

x_strut = 15
x_eng = 7
d_fus = 6.8
b = 25

E_wing = 69 * (10 ** 9)
I_wing = 0.005

l_max = 42580
l_min = 8011
l_slope = (l_min - l_max) / b

w_max = 986
w_min = 95
w_slope = (w_min - w_max) / b

w_eng = 73144

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
    
    length = 0.75 * b
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
        
        deflections.append(defl)
        
    return x_range, deflections
    
plt.plot(x_range, deflections)
plt.show()
print(defl * 1000)