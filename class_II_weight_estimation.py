from constants_and_conversions import g_0, lbs_to_kg, ft_to_m, kts_to_ms, hp_to_N, psf_to_nm2, nm_to_km
import numpy as np


class Class_II:
    
    def __init__(self, W_TO, b, S, n_eng, m_eng, V_D, T_TO):
        self.W_TO = ((W_TO / g_0) / lbs_to_kg)
        self.b = b / ft_to_m
        self.S = S / (ft_to_m ** 2)
        self.V_D = V_D / kts_to_ms
        self.T_TO = ((T_TO / g_0) / lbs_to_kg)
        self.n_eng = n_eng
        self.m_eng = m_eng / lbs_to_kg

    def wing_weight(self, w_f, semi_chord_sweep, n_ult, t_max, choice):
        # w_f is the fuel weight
        # semi_chord_sweep is the semi_chord sweep angle of the wing
        # n_ult is the ultimate load factor
        # t_max is the maximum thickness at the wing root
        
#        print("the wmzf" + str((w_to - w_f)))
#        print("span" + str(b))
#        print((t_max))
#        print(semi_chord_sweep)
        
        w_mzf = self.W_TO - (w_f / lbs_to_kg) / g_0
        t_max = t_max / ft_to_m
        s_angle = np.cos(semi_chord_sweep)
    
        w_weight = 0.0017 * w_mzf * ((self.b / s_angle) ** 0.75) * (1 + np.sqrt((6.3 * s_angle) / self.b)) * (n_ult ** 0.55) * (
                ((self.b * self.S) / (t_max * w_mzf * s_angle)) ** 0.30)
    
        changes = [0.02, -0.05, -0.10, -0.05, -0.30, 0.02]
        # first is to include spoilers and speed brakes
        # second is with 2 wing mounted engines
        # third is with 4 wing mounted engines
        # fourth is if landing gear is not wing mounted
        # fifth is for strutted wings
        # sixth is for fowler flaps
    
        for i in range(len(choice)):
            if choice[i] == 1:
                w_weight = w_weight + w_weight * changes[i]
#                print(w_weight)
    
        return w_weight
    
    
    def _tail_weight(self, k, s, semi_chord_sweep):
        # variables are explained under empennage_weight function
        return k * s * ((3.81 * (s ** 0.2) * self.V_D) / (1000 * (np.cos(semi_chord_sweep) ** 0.5)) - 0.287)
    
    
    def empennage_weight(self, choice, surface, sweep, z_h, span_v):
        # surface is an array of the horizontal and vertical tail surface area
        # v_d is the design dive speed
        # sweep is an array list of the semi_chord sweep angle of the horizontal and vertical tail
        # z_h is the distance between root of vertical tail and start of horizontal tail on vertical tail
        # span_v is the span of the vertical tail

        surface = surface / (ft_to_m ** 2)
        z_h = z_h / ft_to_m
        span_v = span_v / ft_to_m
    
        k_h = 1
        # choose 1 if you have variable incidence stabilizers
        if choice[0] == 1:
            k_h = 1.1
    
        weight_h_tail = self._tail_weight(k_h, surface[0], sweep[0])
    
        k_v = 1
        # choose 1 if the horizontal tails are fin mounted
        if choice[1] == 1:
            k_v = 1 + 0.15 * ((surface[0] * z_h) / (surface[1] * span_v))
    
        weight_v_tail = self._tail_weight(k_v, surface[1], sweep[1])
    
        return weight_h_tail + weight_v_tail
    
    
    def fuselage_weight(self, choice, l_h, w_f, h_f, s_fgs):
        # v_d is the design dive speed
        # l_h is the tailarm from c/4 wing to c/4 tail in ft^2
        # w_f is maximum width of the fuselage
        # f_h is maximum height of the fuselage
        # s_fgs = fuselage fross shell area in ft^2
        
#        print(s_fgs)
#        print(self.V_D)
#        print(w_f)
#        print(h_f)
#        print(l_h)

        w_f = w_f / ft_to_m
        h_f = h_f / ft_to_m
        l_h = l_h / ft_to_m
        s_fgs = s_fgs / (ft_to_m ** 2)
    
        k_f = 1.08
        # choose 1 if the main landing gear is attached to the fuselage
        if choice == 1:
            k_f = k_f * 1.07
    
        return 0.021 * k_f * (((self.V_D * l_h) / (w_f + h_f)) ** 0.5) * (s_fgs ** 1.2)
    
    
    def nacelle_weight(self, choice):
        # t_to is the take-off thrust
        # note that the weight contains ALL nacelles
        # choose 1 if a turbojet or low bypass ratio turbofan is used
        
        if choice == 1:
            return 0.055 * self.T_TO 
        else:
            return 0.065 * self.T_TO
    
    
    def landing_gear_weight(self):
        # w_to is the take-off weight of the aircraft
        
        return 62.21 * ((self.W_TO / 1000) ** 0.84)
    
    
    def engine_weight(self):
        return self.n_eng * self.m_eng
    
    
    @staticmethod
    def induction_weight(l_d, n_inl, a_inl, choice):
        # l_d is the duct length
        # n_inl is the number of inlets
        # a_inl is the cross sectional area of an inlet
        
        if choice[0] == 1:
            return 0
        else:
            l_d = l_d / ft_to_m
            a_inl = a_inl / (ft_to_m ** 2)
            k_d = 1.
    
            # choose 1 if the ducts have a flat cross section
            if choice[1] == 1:
                k_d = 1.33
    
            return 11.45 * ((l_d * n_inl * (a_inl ** 0.5) * k_d) ** 0.7331)
    
    
    @staticmethod
    def propeller_weight(choice, n_p, d_p, p_to, n_bl):
        # n_p is the number of propellers
        # n_bl is the number of blades per propeller
        # p_to is the required take-off power
        # d_p is the propeller diameter
        k_prop = 0.108
        d_p = d_p / ft_to_m
        p_to = p_to / hp_to_N
    
        # choose 1 if piston engines are used
        if choice == 1:
            k_prop = 0.144
    
        return k_prop * (n_p ** 0.218) * ((d_p * p_to * (n_bl ** 0.5)) ** 0.782)
    
    
    def fuel_system_weight(self, n_t, w_f, choice):
        # n_eng is the number of engines
        # n_t is the number of separate fuel tanks
        # k_fsp is for aviation gasoline
        # w_f is the fuel weight
    
        k_fsp = 5.87
        w_f = (w_f / lbs_to_kg) / g_0
    
        # choose 1 if you hae non self-sealing bladder tanks
        if choice == 1:
            return 1.6 * ((w_f / k_fsp) ** 0.727)
        else:
            return 80 * (self.n_eng + n_t - 1) + 15 * (n_t ** 0.5) * ((w_f / k_fsp) ** 0.333)
    
    
    def calc_w_ec(self, l_f, choice):
        # n_eng is the number of engines
        # l_f is the fuselage length
        # b is the span of the main wing
        # the first choice is regarding the type of engine controls: 1 for fuselage mounted jet, 2 for wing mounted jet,
        # 3 for wing mounted turboprops and 4 for wing mounted piston engines
        # the second choice should be 1 if an afterburner is present
    
        l_f = l_f / ft_to_m

        if choice[0] == 1:
            k_ec = 0.686
            if choice[1] == 1:
                k_ec = 1.080
    
            return k_ec * ((l_f * self.n_eng) ** 0.792)
    
        elif choice[0] == 2:
            return 88.46 * ((((l_f + self.b) * self.n_eng) / 100) ** 0.294)
    
        elif choice[0] == 3:
            return 56.84 * ((((l_f + self.b) * self.n_eng) / 100) ** 0.514)
    
        else:
            return 60.27 * ((((l_f + self.b) * self.n_eng) / 100) ** 0.724)
    
    
    def calc_w_ess(self, choice):
        # w_e is the weight per engine
        # n_eng is the number of engines
        # choice is the type of starting system, 1 for one or two jet engines with pneumatic starting system
        # 2 for four jet engines with pneumatic starting systems, 3 for jet engines using electric starting systems
        # 4 for turboprops with pneumatics, 5 for piston engines using electric systems
    
        if choice == 1:
            return 9.33 * (self.m_eng / 1000) ** 1.078
        elif choice == 2:
            return 49.19 * (self.m_eng / 1000) ** 0.541
        elif choice == 3:
            return 38.93 * (self.m_eng / 1000) ** 0.918
    
        elif choice == 4:
            return 12.05 * (self.m_eng / 1000) ** 1.458
    
        elif choice == 5:
            return 50.38 * (self.m_eng / 1000) ** 0.459
    
    
    def calc_w_pc(self, n_bl, n_p, d_p, p_to, choice):
        # choice depends on the typ of engines, 1 is for jets, 2 for turboprops and 3 for piston engines
        # input parameters are as defined previously
        d_p = d_p / ft_to_m
        p_to = p_to / hp_to_N
    
        if choice == 1:
            return 0
        elif choice == 2:
            return 0.322 * (n_bl ** 0.589) * (((n_p * d_p * p_to / self.n_eng) / 1000) ** 1.178)
        else:
            return 4.552 * (n_bl ** 0.379) * (((n_p * d_p * p_to / self.n_eng) / 1000) ** 0.759)
    
    
    def calc_w_osc(self, choice):
        # choice depends on engines, 1 for jet, 2 for turboprop, 3 for radial piston engines, 4 for horizontally opposed
        # piston engines
        # w_e is the weight per engine
        # n_eng is the number of engines
    
        if choice == 1:
            return 0
        elif choice == 2:
            return 0.07 * m_eng
        elif choice == 3:
            return 0.08 * m_eng
        else:
            return 0.03 * m_eng
    
    
    def calc_w_fc(self, q_d):
        # q_d is the dynamic pressure in dive conditions
        # w_to is the take-off weight
    
        q_d = q_d / psf_to_nm2
    
        return 56.01 * (((self.W_TO * q_d) / 100000) ** 0.576)
    
    
    def calc_w_hps_els(self, W_E, choice, v_pax):
        # choice depends on engines, 1 is propellers
        # w_to is the take-off weight
        # v_pax is the passenger cabin volume
    
        v_pax = v_pax / (ft_to_m ** 3)
    
        if choice == 1:
            return 0.325 * (W_E ** 0.8)
        else:
            w_hyd = 0.009 * self.W_TO
            w_els = 10.8 * (v_pax ** 0.7) * (1 - 0.018 * (v_pax ** 0.35))
            return w_hyd + w_els
    
    
    @staticmethod
    def calc_w_instr(w_e, r_max):
        # w_e is the empty weight of the aircraft
        # r is the maximum range of the aircraft
    
        w_e = (w_e / lbs_to_kg) / g_0
        r_max = r_max / nm_to_km
    
        return 0.575 * (w_e ** 0.556) * (r_max ** 0.25)
    
    
    @staticmethod
    def calc_w_api(v_pax, n_crew, n_pax):
        # v_pax is the passenger cabin volume
        v_pax = v_pax / (ft_to_m ** 3)
    
        return 469 * (((v_pax * (n_crew + n_pax)) / 10000) ** 0.419)
    
    
    @staticmethod
    def calc_w_ox(n_crew, n_pax):
        return 7 * ((n_crew + n_pax) ** 0.702)
    
    
    def calc_w_apu(self):
        # w_to is the take-off weight
        return 0.0085 * self.W_TO
    
    
    def calc_w_fur(self, w_f):
        # if more detail is required, look in Torenbeek p292
        # w_to is the take-off weight
        # w_f is the fuel weight of the aircraft
    
        w_f = (w_f / lbs_to_kg) / g_0
    
        return 0.211 * ((self.W_TO - w_f) ** 0.91)
    
    
    @staticmethod
    def calc_w_bc(s_ff):
        # cargo containers require more investigatoin, roskam 110
        # s_ff is the freight floor area
    
        s_ff = s_ff / (ft_to_m ** 2)
    
        return 3 * s_ff
