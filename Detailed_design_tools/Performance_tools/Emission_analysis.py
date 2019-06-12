# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 17:11:54 2019

@author: nikki
"""

import numpy as np
import matplotlib.pytplot as plt

#----------------------------------INPUTS--------------------------------------
#SAR = cruise fuel consumption in kg/km-pax
#N_eng = number of engines

#-----------------------------------STANDARD INPUTS---------------------------
"""Time in mode [s]"""
TIM_TO = 0.7*60. 
TIM_climb = 2.2*60.
TIM_app = 4.*60.
TIM_idle = 26.*60.  #Includes taxiing and holding position

"""Thrust settings as fractions"""
TS_TO = 1.
TS_climb = 0.85
TS_app = 0.3
TS_idle = 0.07
TS_cruise = 0.6

#-----------------------------------DEFINITIONS--------------------------------