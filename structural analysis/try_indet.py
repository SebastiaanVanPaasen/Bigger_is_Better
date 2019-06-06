import numpy as np

x_strut = 15
x_eng = 7
d_fus = 6.8
b = 25

E_wing = 69 * (10 ** 9)
I_wing = 0.005
	
gamma = np.tan(d_fus / x_strut)

A_strut = 0.25 * np.pi * ((4 / 100) ** 2)
L_strut = np.sqrt(d_fus ** 2 + x_strut ** 2)
E_strut = 69 * (10 ** 9)

l_max = 42580
l_min = 8011
l_slope = (l_min - l_max) / b

w_max = 986
w_min = 95
w_slope = (w_min - w_max) / b

w_eng = 73144

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
shear_strut = lift_aft_cst + lift_aft_diff + weight_aft_cst + weight_aft_diff

V_root = lift + weight + shear_strut - strut_force

mom_aft_strut = -1 * (mom_aft_cst + mom_aft_diff)
mom_shear_aft_strut = -1 * shear_strut * x_strut
mom_eng = x_eng * w_eng
mom_lift = -lift_front_cst * (x_strut ** 2) / 2 - 0.5 * lift_front_diff * (x_strut ** 2) / 3
mom_weight = -1 * weight_front_cst * (x_strut ** 2) / 2 + 0.5 * weight_front_diff * (x_strut ** 2) / 3
mom_strut_force = strut_force * x_strut

M_root = -1 * (mom_aft_strut + mom_shear_aft_strut + mom_eng + mom_lift + mom_weight + mom_strut_force)

print(M_root)
print(lift)
print("Sturt force in N: ", strut_force)
print("Wing deflection in mm: ", d_wing * 1000)
print("Strut deflection in mm: ", d_strut * 1000)
print("Strut stress in mpa: ", (strut_force / A_strut) / (10 ** 6))

