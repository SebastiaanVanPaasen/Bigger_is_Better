import numpy as np
import constants_and_conversions as cc

from scipy import interpolate
from class_I.lift_distr import lift_distribution, get_correct_data


class Section:
    
    
    def __init__(self, C_root, C_tip, b, x, dx, tc, n, cd_0):
        self.cr = C_root
        self.ct = C_tip
        self.b = b
        self.x = x
        self.width = dx
        self.tc = tc
        self.n = n
        self.cd_0 = cd_0
        self.chord = self.calc_chord()
        self.forces = np.array([])
        

    def calc_chord(self):
        return self.cr - ((self.cr - self.ct) / (self.b / 2)) * self.x
    
    
    def calc_area(self):
        self.area = self.width * self.chord
        
    
    def calc_volume(self):
        self.volume = self.width * self.chord * (self.tc * self.chord)
        
    
    def calc_sc(self, le_sweep):
        self.sc =  0.4 * self.chord + np.tan(le_sweep) * self.x
        
        
    def run_geometrics(self, le_sweep):
        self.calc_chord()
        self.calc_area()
        self.calc_volume()
        self.calc_sc(le_sweep)


    def _calc_coefficients(self, curve_cl, curve_cd):
        cl = curve_cl(self.x)
        
#        print(self.cd_0)
#        print(curve_cd(self.x))
        
        cd = self.cd_0 + curve_cd(self.x)
        
        return cl, cd
    
    
    def calc_forces(self, w_spec, cl_curve, cd_curve, rho, V):
        cl, cd = self._calc_coefficients(cl_curve, cd_curve)
        
#        print(cd)
        
        lift= 1.5 * self.n * (0.5 * rho * (V ** 2) * self.area * cl)
        drag = -1 * 1.5 * self.n * (0.5 * rho * (V **2) * self.area * cd)
        weight = -1 * self.volume * w_spec
        
        self.forces = np.append(self.forces, [drag, lift, weight])

        
    def calc_fuel_weight(self, w_spec_f, fuel_weight):
        fuel_stored = self.volume * w_spec_f

        if fuel_weight > fuel_stored:            
            fuel_weight = -1 * fuel_stored * cc.g_0
            
        elif fuel_weight > 0:
            fuel_stored = fuel_weight
            fuel_weight = -1 * fuel_stored * cc.g_0
            
        else:
            fuel_stored = 0
            fuel_weight = 0
        
        self.forces = np.append(self.forces, fuel_weight)
        fuel_weight = fuel_weight - fuel_stored
        
        return fuel_weight
        
        
    def calc_strut_force(self, loc, F_strut):
        if self.x == loc:
            strut_force = F_strut
        else:
            strut_force = 0
            
        self.forces = np.append(self.forces, strut_force)
        
        
    def calc_engine_char(self, locations, tot_thrust, w_engine, n_eng):
        
        for i in range(len(locations)):
            if self.x == locations[i]:
                thrust = tot_thrust / n_eng
                w_eng = -1 * w_engine
            else:
                thrust = 0
                w_eng = 0
                
        self.forces = np.append(self.forces, [w_eng, thrust])
        
    
    def _calc_torque_positions(self, positions, le_sweep):
        locations = np.zeros((1, len(positions) + 1))

        for i in range(len(locations[0]) - 1):
            locations[0][i] = positions[i] * self.chord + np.tan(le_sweep) * self.x
        
        locations[0][-1] = positions[-1] * self.tc * self.chord
        
#        print(locations)
        
        return locations
        
    
    def calc_torques(self, positions, le_sweep, sc_sweep, pos_thrust):
        
        positions = self._calc_torque_positions(positions, le_sweep)
        torques = np.zeros((1, len(positions[0])))
        
#        print(self.forces)
#        print(torques)
#        print(positions)
        
        for i in range(1, len(positions[0])):
            torques[0][i - 1] = -self.forces[i] * (positions[0][i - 1] - self.sc) * np.cos(sc_sweep)
            
#            print(-self.forces[i] * (positions[0][i - 1] - self.sc) * np.cos(sc_sweep))
            
        torques[0][-1] = self.forces[-1] * (positions[0][-1] - pos_thrust)
        
        T = np.sum(torques)
        
        return T


    def calc_tot_force(self):
        F_y = np.sum(self.forces[1:-1])
        F_z = self.forces[0] + self.forces[-1]

        return F_y, F_z
    

    def calc_moments(self, x_0, positions):
        M_z = np.zeros((1, len(positions) - 2))
        
        for i in range(1, len(positions) - 1):
            M_z[0][i - 1] = self.forces[i] * (self.x - self.width / 2 - x_0)
            
        M_z = np.sum(M_z)
        M_y = -1 * (self.forces[0] + self.forces[-1]) * (self.x - self.width / 2 - x_0)
        
        return M_z, M_y
        
        
class Helpers:
    
    
    def calc_polars(CL):
        output_avl = lift_distribution(CL)
        x_pos, cl, cdi = get_correct_data(output_avl)

        # Lift Code
        PolyFitCurveCl = interpolate.interp1d(x_pos, cl, kind="cubic", fill_value="extrapolate")

        # Drag Code
        PolyFitCurveidrag = interpolate.interp1d(x_pos, cdi, kind='cubic', fill_value='extrapolate')
        
        return PolyFitCurveCl, PolyFitCurveidrag
    
    
    def input_CL(S, V, rho, W):
        return W / (0.5 * rho * (V ** 2) * S)
    
    
#    def S_cross_section(chord, tc):
#        
#        return chord * (tc * chord)
    
    
        
        
        
       