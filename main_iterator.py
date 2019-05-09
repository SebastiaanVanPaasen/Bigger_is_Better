from inputs import *
from class_I_weight_estimation import class_I

weights = class_I(A, Oswald, S_ratio, C_fe, cruise_range, reserve_range, V_cruise, c_j_cruise)
# print("The take-off mass equals " + str(round(weights[0], 2)))
# print("The fuel mass equals " + str(round(weights[3], 2)))
# print("The payload mass equals " + str(round(weights[2], 2)))
# print("The empty mass equals " + str(round(weights[1], 2)))

design_point = np.array([0.3, 6000])

thrust, surface = weights[0]*design_point[0], weights[0]/design_point[1]





