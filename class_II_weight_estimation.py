from constants_and_conversions import *
import numpy as np


def wing_weight(w_to, w_f, b, semi_chord_sweep, n_ult, s, t_max, choice):
    # w_to is the take-off weight
    # w_f is the fuel weight
    # semi_chord_sweep is the semi_chord sweep angle of the wing
    # n_ult is the ultimate load factor
    # s is the surface area of the wing
    # t_max is the maximum thickness at the wing root
    # print("the wmzf" + str((w_to - w_f)))
    # print("span" + str(b))
    # print((t_max))
    # print(semi_chord_sweep)
    # print(s)
    w_mzf = ((w_to - w_f) / lbs_to_kg) / g_0
    b = b / ft_to_m
    t_max = t_max / ft_to_m
    s_angle = np.cos(semi_chord_sweep)
    s = s / (ft_to_m ** 2)

    w_weight = 0.0017 * w_mzf * ((b / s_angle) ** 0.75) * (1 + np.sqrt(6.3 * s_angle / b)) * (n_ult ** 0.55) * (
            ((b * s) / (t_max * w_mzf * s_angle)) ** 0.30)

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

    return w_weight


def _tail_weight(k, s, v_d, semi_chord_sweep):
    # variables are explained under empennage_weight function
    return k * s * ((3.81 * (s ** 0.2) * v_d) / (1000 * np.cos(semi_chord_sweep) ** 0.5) - 0.287)


def empennage_weight(choice, surface, v_d, sweep, z_h, span_v):
    # surface is an array of the horizontal and vertical tail surface area
    # v_d is the design dive speed
    # sweep is an array list of the semi_chord sweep angle of the horizontal and vertical tail
    # z_h is the distance between root of vertical tail and start of horizontal tail on vertical tail
    # span_v is the span of the vertical tail
    v_d = v_d / kts_to_ms
    surface = surface / (ft_to_m ** 2)
    z_h = z_h / ft_to_m
    span_v = span_v / ft_to_m

    k_h = 1
    # choose 1 if you have variable incidence stabilizers
    if choice[0] == 1:
        k_h = 1.1

    weight_h_tail = _tail_weight(k_h, surface[0], v_d, sweep[0])

    k_v = 1
    # choose 1 if the horizontal tails are fin mounted
    if choice[1] == 1:
        k_v = 1 + 0.15 * ((surface[0] * z_h) / (surface[1] * span_v))

    weight_v_tail = _tail_weight(k_v, surface[1], v_d, sweep[1])

    return weight_h_tail + weight_v_tail


def fuselage_weight(choice, v_d, l_h, w_f, h_f, s_fgs):
    # v_d is the design dive speed
    # l_h is the tailarm from c/4 wing to c/4 tail in ft^2
    # w_f is maximum width of the fuselage
    # f_h is maximum height of the fuselage
    # s_fgs = fuselage fross shell area in ft^2
    # print(s_fgs)
    # print(v_d)
    # print(w_f)
    # print(h_f)
    # print(l_h)
    v_d = v_d / kts_to_ms
    w_f = w_f / ft_to_m
    h_f = h_f / ft_to_m
    l_h = l_h / ft_to_m
    s_fgs = s_fgs / (ft_to_m ** 2)
    # print(s_fgs)
    # print(v_d)
    # print(w_f)
    # print(h_f)
    # print(l_h)

    k_f = 1.08
    # choose 1 if the main landing gear is attached to the fuselage
    if choice == 1:
        k_f = k_f * 1.07

    return 0.021 * k_f * (((v_d * l_h) / (w_f + h_f)) ** 0.5) * (s_fgs ** 1.1)


def nacelle_weight(t_to, choice):
    # t_to is the take-off thrust
    # note that the weight contains ALL nacelles
    # choose 1 if a turbojet or low bypass ratio turbofan is used
    if choice == 1:
        return 0.055 * ((t_to / lbs_to_kg) / g_0)
    else:
        return 0.065 * ((t_to / lbs_to_kg) / g_0)


def landing_gear_weight(w_to):
    # w_to is the take-off weight of the aircraft
    return 62.21 * ((((w_to / lbs_to_kg) / g_0) / 1000) ** 0.84)


def engine_weight(n_e, w_e):
    # w_e is the weight of one engine
    w_e = (w_e / lbs_to_kg)
    return n_e * w_e


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


def fuel_system_weight(n_e, n_t, w_f, choice):
    # n_e is the number of engines
    # n_t is the number of separate fuel tanks
    # k_fsp is for aviation gasoline
    # w_f is the fuel weight

    k_fsp = 5.87
    w_f = (w_f / lbs_to_kg) / g_0

    # choose 1 if you hae non self-sealing bladder tanks
    if choice == 1:
        return 1.6 * ((w_f / k_fsp) ** 0.727)
    else:
        return 80 * (n_e + n_t - 1) + 15 * (n_t ** 0.5) * ((w_f / k_fsp) ** 0.333)


def calc_w_ec(l_f, n_e, b, choice):
    # n_e is the number of engines
    # l_f is the fuselage length
    # b is the span of the main wing
    # the first choice is regarding the type of engine controls: 1 for fuselage mounted jet, 2 for wing mounted jet,
    # 3 for wing mounted turboprops and 4 for wing mounted piston engines
    # the second choice should be 1 if an afterburner is present

    l_f = l_f / ft_to_m
    b = b / ft_to_m

    if choice[0] == 1:
        k_ec = 0.686
        if choice[1] == 1:
            k_ec = 1.080

        return k_ec * ((l_f * n_e) ** 0.792)

    elif choice[0] == 2:
        return 88.46 * (((l_f + b) / n_e / 100) ** 0.294)

    elif choice[0] == 3:
        return 56.84 * (((l_f + b) / n_e / 100) ** 0.514)

    else:
        return 60.27 * (((l_f + b) / n_e / 100) ** 0.724)


def calc_w_ess(w_e, n_e, choice):
    # w_e is the weight per engine
    # n_e is the number of engines
    # choice is the type of starting system, 1 for one or two jet engines with pneumatic starting system
    # 2 for four jet engines with pneumatic starting systems, 3 for jet engines using electric starting systems
    # 4 for turboprops with pneumatics, 5 for piston engines using electric systems
    w_e = ((w_e * n_e) / lbs_to_kg)

    if choice == 1:
        return 9.33 * (w_e / 1000) ** 1.078
    elif choice == 2:
        return 49.19 * (w_e / 1000) ** 0.541
    elif choice == 3:
        return 38.93 * (w_e / 1000) ** 0.918

    elif choice == 4:
        return 12.05 * (w_e / 1000) ** 1.458

    elif choice == 5:
        return 50.38 * (w_e / 1000) ** 0.459


def calc_w_pc(n_bl, n_p, d_p, n_e, p_to, choice):
    # choice depends on the typ of engines, 1 is for jets, 2 for turboprops and 3 for piston engines
    # input parameters are as defined previously
    d_p = d_p / ft_to_m
    p_to = p_to / hp_to_N

    if choice == 1:
        return 0
    elif choice == 2:
        return 0.322 * (n_bl ** 0.589) * (((n_p * d_p * p_to / n_e) / 1000) ** 1.178)
    else:
        return 4.552 * (n_bl ** 0.379) * (((n_p * d_p * p_to / n_e) / 1000) ** 0.759)


def calc_w_osc(choice, w_e, n_e):
    # choice depends on engines, 1 for jet, 2 for turboprop, 3 for radial piston engines, 4 for horizontally opposed
    # piston engines
    # w_e is the weight per engine
    # n_e is the number of engines
    w_e = ((w_e * n_e) / lbs_to_kg)

    if choice == 1:
        return 0
    elif choice == 2:
        return 0.07 * w_e
    elif choice == 3:
        return 0.08 * w_e
    else:
        return 0.03 * w_e


def calc_w_fc(w_to, q_d):
    # q_d is the dynamic pressure in dive conditions
    # w_to is the take-off weight

    w_to = (w_to / lbs_to_kg) / g_0
    q_d = q_d / psf_to_nm2

    return 56.01 * (((w_to * q_d) / 100000) ** 0.576)


def calc_w_hps_els(choice, w_to, v_pax):
    # choice depends on engines, 1 is propellers
    # w_to is the take-off weight
    # v_pax is the passenger cabin volume

    w_to = (w_to / lbs_to_kg) / g_0
    v_pax = v_pax / (ft_to_m ** 3)

    if choice == 1:
        return 0.325 * (w_to ** 0.8)
    else:
        w_hyd = 0.009 * w_to
        w_els = 10.8 * (v_pax ** 0.7) * (1 - 0.018 * (v_pax ** 0.35))
        return w_hyd + w_els


def calc_w_instr(w_e, r_max):
    # w_e is the empty weight of the aircraft
    # r is the maximum range of the aircraft

    w_e = (w_e / lbs_to_kg) / g_0
    r_max = r_max / nm_to_km

    return 0.575 * (w_e ** 0.556) * (r_max ** 0.25)


def calc_w_api(v_pax, n_crew, n_pax):
    # v_pax is the passenger cabin volume
    v_pax = v_pax / (ft_to_m ** 3)

    return 469 * (((v_pax * (n_crew + n_pax)) / 10000) ** 0.419)


def calc_w_ox(n_crew, n_pax):
    return 7 * ((n_crew + n_pax) ** 0.702)


def calc_w_apu(w_to):
    # w_to is the take-off weight
    w_to = (w_to / lbs_to_kg) / g_0

    return 0.0085 * w_to


def calc_w_fur(w_to, w_f):
    # if more detail is required, look in Torenbeek p292
    # w_to is the take-off weight
    # w_f is the fuel weight of the aircraft

    w_to = (w_to / lbs_to_kg) / g_0
    w_f = (w_f / lbs_to_kg) / g_0

    return 0.211 * ((w_to - w_f) ** 0.91)


def calc_w_bc(s_ff):
    # cargo containers require more investigatoin, roskam 110
    # s_ff is the freight floor area

    s_ff = s_ff / (ft_to_m ** 2)

    return 3 * s_ff
