import numpy as np


class Section:
    def __init__(self, C_root, C_tip, b, x, dx, tc, n):
        self.cr = C_root
        self.ct = C_tip
        self.b = b
        self.x = x
        self.width = dx
        self.tc = tc
        self.n = n
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


    def _calc_coefficients(self, curve_cl, curve_cd):
        cl = curve_cl(self.x)
        cd = self.cd_0 + curve_cd(self.x)
        
        return cl, cd
    
    
    def calc_forces(self, w_spec, cl_curve, cd_curve, rho, V):
        cl, cd = self._calc_coefficients(cl_curve, cd_curve)
        
        lift= 1.5 * self.n * (0.5 * rho * (V ** 2) * self.area * cl)
        drag = -1 * 1.5 * self.n * (0.5 * rho * (V **2) * self.area * cd)
        weight = -1 * self.volume * w_spec
        self.forces = np.append(self.forces, [lift, drag, weight])

        
    def calc_fuel_weight(self, ff, w_spec_f, limit):
        dist = limit - self.x
        if dist > self.width:
            fuel_weight = -1 * self.volume * w_spec_f * ff
        elif dist > 0:
            part_volume = (self.volume / self.width ) * dist
            fuel_weight = -1 * part_volume * w_spec_f * ff
        else:
            fuel_weight = 0
        
        self.forces = np.append(self.forces, fuel_weight)
        
        
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
        
    
    def calc_positions(self, positions, le_sweep):
        locations = np.zeros((1, len(positions) + 1))
        
        for i in range(len(positions) - 1):
            locations[i] = positions[i] * self.chord + np.tan(le_sweep) * self.x
            
        locations[-1] = positions[-1] * self.tc * self.chord
        
    
    def calc_torques(self, positions, sc_sweep, pos_thrust):
        torques = np.zeros((1, len(positions)))
        
        for i in range(len(positions)):
            torques[i] = -self.forces[i] * (positions[i] - self.sc) * np.cos(sc_sweep)
            
        torques[-1] = self.forces[-1] * (positions[-1] - pos_thrust)
        
        T = np.sum(torques)
        return T


    def calc_tot_force(self):
        v_force = np.sum(self.forces[:-1])
        h_force = self.forces[-1]

        return v_force, h_force
    

    def calc_moments(self, x_0, positions):
        M_z = np.zeros((1, len(positions) - 1))
        
        for i in range(len(positions) - 1):
            M_z[i] = self.forces[i] * (self.x - x_0)
            
        M_z = np.sum(M_z)
        M_y = -1 * self.forces[-1] * (self.x - x_0)
        
        return M_z, M_y
        
        
        
        
        
        
        
       