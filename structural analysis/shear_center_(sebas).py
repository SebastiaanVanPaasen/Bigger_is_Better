import numpy as np
import matplotlib.pyplot as plt
import centroid_wing as cw
import airfoil_geometry as ag
import practised_ultimate as pu
import Airfoil_inertia as ai

print("------- Starting on shear calculation -------")

boom_area = pu.boom_area_all[0]
I_zz, I_yy, I_yz = pu.I_zz_sections, pu.I_yy_wing, pu.I_yz_wing

#print("The boom area", boom_area)

data_z_all_sec, data_y_upper_all_sec, data_y_lower_all_sec = ag.airfoil_geometry(pu.N, pu.b, pu.calc_chord, pu.X_root)
airfoil_area, z_c_airfoil, y_c_airfoil = cw.get_skin_centroid(pu.N, pu.b, pu.calc_chord, pu.dx)
z_centroid_all_sec, y_centroid_all_sec, y_loc_spar_up, y_loc_spar_low, y_loc_stiff_up, y_loc_stiff_low, y_vertical_spar, z_loc_stiff_up, spar_loc_sec, z_loc_stiff_low, spar_areas_verti = cw.wing_centroid(boom_area, cw.spar_areas_hori, cw.t_spar_v, z_c_airfoil, y_c_airfoil, cw.n_stiff_up, cw.n_stiff_low, pu.N, pu.b, pu.calc_chord, pu.X_root, pu.dx)
s_all_sec, ds_all_sec = ai.s_airfoil(pu.N, pu.b, pu.calc_chord, pu.X_root)


boom_loc_y = np.zeros((len(y_loc_spar_up), len(y_loc_stiff_up[0]) * 2 + 4))
boom_loc_z = np.zeros((len(y_loc_spar_up), len(y_loc_stiff_up[0]) * 2 + 4))

for i in range(len(y_loc_spar_up)):
    
    for j in range(len(y_loc_stiff_up[0])):
        boom_loc_y[i][j] = -(y_loc_stiff_up[i][j] - y_centroid_all_sec[i])
        boom_loc_z[i][j] = z_loc_stiff_up[i][j] - z_centroid_all_sec[i]
        
        
    boom_loc_y[i][j + 1] = -(y_loc_spar_up[i][1] - y_centroid_all_sec[i])
    boom_loc_z[i][j + 1] = 0.55 * data_z_all_sec[i][-1] - z_centroid_all_sec[i]
    
    
    boom_loc_y[i][j + 2] = -(y_loc_spar_low[i][1]  - y_centroid_all_sec[i])
    boom_loc_z[i][j + 2] = 0.55 * data_z_all_sec[i][-1] - z_centroid_all_sec[i]
        
    rev_y = y_loc_stiff_low[i][::-1]
    rev_z = z_loc_stiff_low[i][::-1]
    
    for k in range(len(y_loc_stiff_low[0])):
        boom_loc_y[i][j + k + 3] = -(rev_y[k] - y_centroid_all_sec[i])
        boom_loc_z[i][j + k + 3] = rev_z[k] - z_centroid_all_sec[i]
        
        
    boom_loc_y[i][j + k + 4] = -(y_loc_spar_low[i][0] - y_centroid_all_sec[i])
    boom_loc_z[i][j + k + 4] = 0.15 * data_z_all_sec[i][-1] - z_centroid_all_sec[i]

    
    boom_loc_y[i][j + k + 5] = -(y_loc_spar_up[i][0] - y_centroid_all_sec[i])
    boom_loc_z[i][j + k + 5] = 0.15 * data_z_all_sec[i][-1] - z_centroid_all_sec[i]
   

d = []
for i in range(len(boom_loc_y)):
    distances = []
    start_idx = 0
    
    
    start_idx = np.argmin(abs(boom_loc_z[i][0] - data_z_all_sec[i]))
    for j in range(1, len(y_loc_stiff_up[0])):
        
        idx = np.argmin(abs(boom_loc_z[i][j] - data_z_all_sec[i]))
        distance = ds_all_sec[i][start_idx:idx + 1]
        
#        print(distance, start_idx, idx)
        distances.append(sum(distance))
        
        start_idx = idx
        
    idx = np.argmin(abs(boom_loc_z[i][j + 1] - data_z_all_sec[i]))
    distance = ds_all_sec[i][start_idx:idx + 1]
    distances.append(sum(distance))
    
    distances.append(abs(boom_loc_y[i][j + 2] - boom_loc_y[i][j + 1]))
    
    
    start_idx = np.argmin(abs(boom_loc_z[i][j + 2] - data_z_all_sec[i]))
    for k in range(1, len(y_loc_stiff_up[0])):

        idx = np.argmin(abs(boom_loc_z[i][j + k + 2] - data_z_all_sec[i]))
        distance = ds_all_sec[i][idx - 1:start_idx]
        
#        print(distance, start_idx, idx)
        distances.append(sum(distance))
        
        start_idx = idx
        
    idx = np.argmin(abs(boom_loc_z[i][j + k + 3] - data_z_all_sec[i]))
    distance = ds_all_sec[i][idx - 1:start_idx]
    distances.append(sum(distance))
    
    distances.append(abs(boom_loc_y[i][j + k + 3] - boom_loc_y[i][j + k + 4]))
    
    d.append(distances)
    
    
#print(distances)
        
coef_z = I_yz[0] / (I_zz[0] * I_yy[0] - I_yz[0] ** 2)
coef_y = -I_yy[0] / (I_zz[0] * I_yy[0] - I_yz[0] ** 2)   

for i in range(1):
    
    A_spar_front_up = (cw.spar_areas_hori[0] + cw.t_spar_v * (y_loc_spar_up[i][0] - y_loc_spar_low[i][0]) / 6) * (2 + y_loc_spar_low[i][0] / y_loc_spar_up[i][0])
    A_spar_front_low = (cw.spar_areas_hori[0] + cw.t_spar_v * (y_loc_spar_up[i][0] - y_loc_spar_low[i][0]) / 6) * (2 + y_loc_spar_up[i][0] / y_loc_spar_low[i][0])
    front_spar = [A_spar_front_low, A_spar_front_up]
    
    A_spar_back_up = (cw.spar_areas_hori[1] + cw.t_spar_v * (y_loc_spar_up[i][1] - y_loc_spar_low[i][1]) / 6) * (2 + y_loc_spar_low[i][1] / y_loc_spar_up[i][1])
    A_spar_back_low = (cw.spar_areas_hori[1] + cw.t_spar_v * (y_loc_spar_up[i][1] - y_loc_spar_low[i][1]) / 6) * (2 + y_loc_spar_up[i][1] / y_loc_spar_low[i][1])
    back_spar = [A_spar_back_up, A_spar_back_low]
    
    boom_areas = np.append(len(y_loc_stiff_up[0]) * [boom_area], back_spar)
    boom_areas = np.append(boom_areas, len(y_loc_stiff_up[0]) * [boom_area])
    boom_areas = np.append(boom_areas, front_spar)
    
    diff_q = np.zeros(len(boom_loc_y[0]))
    base_q = np.zeros(len(boom_loc_y[0])+1)
    for j in range(len(boom_loc_y[0])):

        delta_q_z = coef_z * boom_loc_z[i][j] * boom_areas[j]
        delta_q_y = coef_y * boom_loc_y[i][j] * boom_areas[j]
        
        print(delta_q_z, delta_q_y)
        diff_q[j] = delta_q_z + delta_q_y
    
    
#    previous = diff_q[0]
    for j in range(len(diff_q)+1):
        base_q[j] = sum(diff_q[:j])
        
#        previous = diff_q[j]

print(diff_q)
print(base_q)
#for i in range(len(boom_locations)):
    
#plt.plot(data_z_all_sec[0], data_y_upper_all_sec[0])
#plt.plot(data_z_all_sec[0], data_y_lower_all_sec[0])
#plt.axis("scaled")
#plt.show()

def q_b(y_loc, z_loc, boom_area, moi):
    """
    Calculates the base shear flow due to a dummy force of 1 N.

    :param boom_locations: List of float lists, containing all the locations of the booms around the cross-section.
    :param area: List of floats, containing the areas of all the booms around the cross-section.
    :param moi: Float, representing the moment of inertia around the z-axis.
    :return: List of floats, containing the base shear flows around the cross-section.
    """
    
    diff_q = np.zeros(len(y_loc))
    base_q = np.zeros(len(y_loc))

    I_zz = moi[0]
    I_yy = moi[1]
    I_zy = moi[2]

    coef_z = -1 / (I_yy - I_zy ** 2)
    coef_y = -1 / (I_zz - I_zy ** 2)
    
    for i in range(len(y_loc)):
        delta_q_z = coef_z * z_loc[i] * boom_area
        delta_q_y = coef_y * y_loc[i] * boom_area
        
        diff_q[i] = delta_q_z + delta_q_y


    previous = diff_q[0]
    
    
    for i in range(1, len(diff_q)):
        base_q[i] = previous + diff_q[i]
        
        previous = diff_q[i]

    # print(base_q)
    # print(len(base_q))
    
    return base_q


def q_s0(base_q, t, s):
    """
    Calculates the constant shear flow equations leaving them as unknowns in the equations.

    :param base_q: List of floats, containing the base shear flows between all the booms in one cell.
    :param t: List of floats, containing the thickness between all the booms in one cell.
    :param s: List of floats, containing the distances between the booms in one cell.
    :return: values for the two unknown shear flows and the numerical value due to the base shear flow.
    """
    flow = np.multiply(base_q, s)
    flow = np.divide(flow, t)

    num = np.sum(flow)
    x = np.sum(np.divide(s, t))
    y = -1 * (s[-1] / t[-1])

    return x, y, num


#def cell_separation(skin, base_q, s, n_top, n_bottom, d_spar, spar, h_a, c_a, booms):
#    """
#    Creates the cell parameters for both cells.
#
#    :param skin: Float, representing the thickness of the skin.
#    :param base_q: List of floats, representing all the base shear flows around the section.
#    :param s: List of floats, representing the circumferential positions of the booms around the section.
#    :param n_top: Integer, boom number of the boom at the top of the spar.
#    :param n_bottom: Integer, boom number of the boom at the bottom of the spar.
#    :param d_spar: Float, circumferential distance to the boom at the top of the spar.
#    :param spar: Float, representing the thickness of the spar.
#    :param h_a: Float, representing the height of the spar.
#    :param c_a: Float, representing the length of the chord.
#    :param booms: List of float lists, containing the locations of all the booms in the yz-coordinate frame.
#    :return: Two lists, containing the base shear flow, circumferential distances between booms, thickness of the
#    distances between the booms, the area of the cell and the coordinates of all the booms in the cell.
#    """
#    base_q_cell_i = np.append(base_q[n_bottom + 1:], base_q[:n_top + 1])
#    base_q_cell_i = np.append(base_q_cell_i, np.array([0]))
#    base_q_cell_ii = np.append(base_q[n_top + 1:n_bottom + 1], np.array([0]))
#
#    bottom_locations = s[n_bottom:]
#    top_locations = s[:n_top + 1]
#
#    s_cell_i_bottom = np.zeros(len(bottom_locations))
#    for i in range(len(bottom_locations) - 1):
#        s_cell_i_bottom[i] = bottom_locations[i + 1] - bottom_locations[i]
#    s_cell_i_bottom[-1] += top_locations[0] * 2
#
#    s_cell_i_top = np.zeros(len(top_locations) - 1)
#    for i in range(len(top_locations) - 1):
#        s_cell_i_top[i] = top_locations[i + 1] - top_locations[i]
#
#    s_cell_i = np.append(np.append(s_cell_i_bottom, s_cell_i_top), np.array([d_spar]))
#
#    s_cell_ii_locations = s[n_top:n_bottom + 1]
#    s_cell_ii = np.zeros(len(s_cell_ii_locations) - 1)
#    for i in range(len(s_cell_ii_locations) - 1):
#        s_cell_ii[i] = s_cell_ii_locations[i + 1] - s_cell_ii_locations[i]
#
#    s_cell_ii = np.append(s_cell_ii, np.array([d_spar]))
#
#    t_cell_i = np.append(np.array((len(s_cell_i) - 1) * [skin]), np.array([spar]))
#    t_cell_ii = np.append(np.array((len(s_cell_ii) - 1) * [skin]), np.array([spar]))
#
#    cell_i_coordinates = booms[n_bottom:] + booms[:n_top + 1]
#    cell_ii_coordinates = booms[n_top:n_bottom + 1]
#
#    area_i = np.pi * ((h_a / 2) ** 2) / 2
#    area_ii = h_a * (c_a - (h_a / 2))
#
#    return [base_q_cell_i, s_cell_i, t_cell_i, cell_i_coordinates, area_i], \
#           [base_q_cell_ii, s_cell_ii, t_cell_ii, cell_ii_coordinates, area_ii]
#
#
#def area_triangle(pos1, pos2, pos3):
#    len1 = np.sqrt((pos1[0] - pos2[0]) ** 2 + (pos1[1] - pos2[1]) ** 2)
#    len2 = np.sqrt((pos2[0] - pos3[0]) ** 2 + (pos2[1] - pos3[1]) ** 2)
#    len3 = np.sqrt((pos1[0] - pos3[0]) ** 2 + (pos1[1] - pos3[1]) ** 2)
#
#    s = (len1 + len2 + len3) / 2
#    area = np.sqrt((s * (s - len1) * (s - len2) * (s - len3)))
#
#    return area
#
#
#def area_circle(distance, h_a):
#    ratio = distance / (h_a * np.pi / 2)
#
#    area = ratio * np.pi * ((h_a / 2) ** 2) / 2
#
#    return area
#
#
#def sum_moments_i(base_flow, constant_flow, cell_area, s, h_a):
#    moment = 0
#
#    for i in range(len(s) - 1):
#        area = area_circle(s[i], h_a)
#        moment += 2 * area * base_flow[i]
#
#    moment += 2 * cell_area * constant_flow
#
#    return moment
#
#
#def sum_moments_ii(base_flow, constant_flow, cell_area, pos, summation_point, c_a):
#    moment = 0
#    flow_areas = []
#
#    for i in range(len(pos) - 1):
#        flow_area = area_triangle(summation_point, pos[i], pos[i + 1])
#        flow_areas.append(flow_area)
#
#    midpoint = len(flow_areas) // 2
#
#    flow_areas[midpoint] += area_triangle(pos[midpoint], pos[midpoint + 1], [c_a, 0])
#
#    for i in range(len(base_flow) - 1):
#        moment += 2 * flow_areas[i] * base_flow[i]
#
#    moment += 2 * constant_flow * cell_area
#
#    return moment
#
#
#def get_shear_flows(b, areas, i_zz, t_skin, d, top, bottom, h_a, t_spar, c_a):
#    q_diff = q_b(b, areas, i_zz)
#
#    cell_i, cell_ii = cell_separation(t_skin, q_diff, d, top, bottom, h_a, t_spar, h_a, c_a, b)
#
#    qs_i = np.array([q_s0(cell_i[0], cell_i[2], cell_i[1])])
#    qs_ii = np.array([q_s0(cell_ii[0], cell_ii[2], cell_ii[1])])
#
#    sys = np.array([[qs_i[0][0], qs_i[0][1]], [qs_ii[0][0], qs_ii[0][1]]])
#    sys_vec = np.array([[-qs_i[0][2]], [-qs_ii[0][2]]])
#
#    values = np.linalg.solve(sys, sys_vec)
#
#    return values, cell_i, cell_ii
#
#
#def summing_moments(constants, cell_i, cell_ii, z_pos):
#    moments_i = sum_moments_i(cell_i[0], constants[0][0], cell_i[4], cell_i[1], height)
#    moments_ii = sum_moments_ii(cell_ii[0], constants[1][0], cell_ii[4], cell_ii[3], [height / 2, 0], chord)
#
#    z_pos_sc = ((moments_i + moments_ii) / -1)
#
#    return z_pos_sc
#
#
## pre shear-flow calculations
#n_booms = 260
#b, boom_distances, boom_areas, top_spar, bottom_spar = get_boom_information(n_booms)
#z_pos, boom_loc = get_cg(n_booms)
#
#izz, iyy = get_moi(n_booms)
## required parameters
#height = 0.173
#chord = 0.484
#skin_t = 1.1 / 1000
#spar_t = 2.5 / 1000
#
#
#def get_shear_center():
#    constant, first_cell, second_cell = get_shear_flows(boom_loc, boom_areas, izz, skin_t, boom_distances,
#                                                        top_spar,
#                                                        bottom_spar, height, spar_t, chord)
#
#    shear_center = summing_moments(constant, first_cell, second_cell, z_pos)
#
#    return shear_center
#
#
#sc_z = get_shear_center()
## print(sc_z)
