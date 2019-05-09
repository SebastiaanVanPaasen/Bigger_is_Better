import numpy as np
# Lift Estimation
# VALUES ARE USED FOR FIRST ESTIMATE USE THE ACTUAL INPUTS!


# Function to find CL_alpha, CL and alpha trim in the case of no wing twist
def Clean_Wing_Lift(aspect_ratio, m_cruise, half_chord_sweep, quarter_chord_sweep, alpha, alpha_0, cl_des, airfoil_eff):
    beta = np.sqrt(1 - m_cruise)  # prandtl_glauert correction

    # DATCOM method for CL_alpha
    CL_alpha = (2 * np.pi * aspect_ratio) / (
            2 + np.sqrt(4 + (aspect_ratio * beta / airfoil_eff) ** 2 * (1 + np.tan(half_chord_sweep) ** 2 / beta ** 2)))

    CL = CL_alpha * (alpha - alpha_0)  # calculate the Lift
    CL_des = cl_des / np.cos(quarter_chord_sweep) ** 2  # CL_des based on Cl_des from airfoil
    alpha_0_L = alpha_0  # In case of no wing twist
    alpha_trim = CL_des / CL_alpha + alpha_0_L  # Find trim angle to obtain required lift

    return np.array([CL_alpha, CL, alpha_trim])


A = 9
M_cruise = 0.7
half_chord_sweep = 0.5
quarter_chord_sweep = 0.6
alpha = 0.5
alpha_0 = 0.3
Cl_des = 1.2
airfoil_eff_factor = 0.95
