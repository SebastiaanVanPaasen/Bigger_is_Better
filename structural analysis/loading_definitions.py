import numpy as np
import constants_and_conversions as cc
import matplotlib.pyplot as plt

from scipy import interpolate
from class_I.lift_distr import lift_distribution, get_correct_data


class Section:
    
    
    def __init__(self, C_root, C_tip, b, x, dx, tc, n, cd_0, pos_eng):
        self.cr = C_root
        self.ct = C_tip
        self.b = b
        self.x = x
        self.width = dx
        self.tc = tc
        self.n = n
        self.cd_0 = cd_0
        self.chord = self.calc_chord()
        self.pos_eng = pos_eng
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
        
        cd = self.cd_0 + curve_cd(self.x) * self.n * 1.5
        
        return cl, cd
    
    
    def calc_forces(self, w_spec, cl_curve, cd_curve, rho, V):
        cl, cd = self._calc_coefficients(cl_curve, cd_curve)
        
#        print(cd)
        
        lift= 1.5 * self.n * (0.5 * rho * (V ** 2) * self.area * cl)
        drag = -1 * (0.5 * rho * (V **2) * self.area * cd)
        weight = -1 * self.volume * w_spec
        
        self.forces = np.append(self.forces, [drag, lift, weight])
        
        return drag

        
    def calc_fuel_weight(self, w_spec_f, fuel_mass):
        fuel_stored = self.volume * w_spec_f
        
#        print("fuel weight", fuel_mass)
#        print("fuel stored", fuel_stored)
#        print("fuel volume", self.volume)
        
        if fuel_mass > fuel_stored:            
            fuel_weight = -1 * fuel_stored * cc.g_0
            
        elif fuel_mass > 0:
            fuel_stored = fuel_mass
            fuel_weight = -1 * fuel_stored * cc.g_0
            
        else:
#            print("at section", self.x)
            fuel_stored = 0
            fuel_weight = 0
        
        self.forces = np.append(self.forces, fuel_weight)
        fuel_mass = fuel_mass - fuel_stored
        
        return fuel_mass
        
        
    def calc_strut_force(self, loc, F_strut):
        if self.x - self.width / 2 <= loc < self.x + self.width / 2:
#            print("strut has been found")
            strut_force = -F_strut
        else:
            strut_force = 0
            
        self.forces = np.append(self.forces, strut_force)
        
        
    def calc_engine_char(self, locations, tot_thrust, w_engine, n_eng):
        
        for i in range(len(locations)):
            if self.x - self.width / 2 <= locations[i] < self.x + self.width / 2:
#                print("found the engine")
                thrust = tot_thrust / n_eng
                w_eng = -1 * w_engine
               
            else:
                thrust = 0
                w_eng = 0

        self.forces = np.append(self.forces, [w_eng, thrust])
        
        
    def _calc_torque_distances(self, positions, le_sweep):
        distances = np.zeros((1, len(positions) + 1))

        for i in range(len(distances[0]) - 1):
            distances[0][i] = positions[i] * self.chord + np.tan(le_sweep) * self.x
        
        distances[0][-1] = positions[-1] * self.tc * self.chord
        
#        print(locations)
        
        return distances
        
    
    def calc_torques(self, positions, le_sweep, sc_sweep, pos_thrust):
        
        distances = self._calc_torque_distances(positions, le_sweep)
        torques = np.zeros((1, len(distances[0])))
        
#        print(self.forces)
#        print(torques)
#        print(positions)
        
        for i in range(1, len(distances[0])):
            torques[0][i - 1] = -self.forces[i] * (distances[0][i - 1] - self.sc) * np.cos(sc_sweep)
            
#            print(-self.forces[i] * (positions[0][i - 1] - self.sc) * np.cos(sc_sweep))
            
        torques[0][-1] = self.forces[-1] * (distances[0][-1] - pos_thrust)
        
        T = np.sum(torques)
        
        return T


    def calc_tot_force(self):
        F_y = np.sum(self.forces[1:-1])
        F_z = self.forces[0] + self.forces[-1]
        
#        print("drag", self.forces[0])
#        print("thrust", self.forces[-1])

        return F_y, F_z
    

    def calc_moments(self):
        M_z = np.zeros((1, len(self.forces) - 2))
        
        for i in range(1, len(self.forces) - 1):
            M_z[0][i - 1] = self.forces[i] * self.x
            
        M_z = np.sum(M_z)
        
        M_y = -1 * self.forces[0] * self.x - 1 * self.pos_eng * self.forces[-1]

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
    
    
    def plotter(discr, T_plot, M_z_plot, M_y_plot, F_y_plot, F_z_plot):
        plt.subplot(2,3,5)
        plt.subplot(2,3,1)
        plt.gca().set_title('Mz distribution')
        plt.plot(discr, M_z_plot)
        plt.xlabel('Position Along Wing Span [$m$]')
        plt.ylabel('Mz $[Nm]$')
        plt.subplot(2,3,2)
        plt.gca().set_title('Fz distribution')
        plt.plot(discr, F_z_plot)
        plt.xlabel('Position Along Wing Span [$m$]')
        plt.ylabel('Fz [$N$]')
        plt.subplot(2,3,3)
        plt.gca().set_title('My distribution')
        plt.plot(discr, M_y_plot)
        plt.xlabel('Position Along Wing Span [$m$]')
        plt.ylabel('My $[Nm]$')
        plt.subplot(2,3,4)
        plt.gca().set_title('Fy distribution')
        plt.plot(discr, F_y_plot)
        plt.xlabel('Position Along Wing Span [$m$]')
        plt.ylabel('Fy [$N$]')
        plt.subplot(2,3,5)
        plt.gca().set_title('T distribution')
        plt.plot(discr, T_plot)
        plt.xlabel('Position Along Wing Span [$m$]')
        plt.ylabel('T $[Nm]$')
        plt.show()
    
    
#    def S_cross_section(chord, tc):
#        
#        return chord * (tc * chord)
    
    
        
        
        
       