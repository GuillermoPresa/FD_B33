#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 15:04:57 2020

@author: frederik
"""

###### Data Processing tasks ########

import math as m
import numpy as np
import matplotlib.pyplot as plt
from cg_pos import cg, payload_list
import ISA_calculator
import Data_reader

#### task 1: First Stationary Measurement Series

# General inputs
gamma = 1.4                                     # [-]
lam = -0.0065                                   # [-]
Ws = 60500                                      # [N] Reduced aircraft weight from assingment
R = 287.0                                       # [J/(mol*K)] Specific gas Constant of Air

# Sea level conditions
p0 = ISA_calculator.pressure0                   # [N/m^2] Pressure at Sea level
g0 = ISA_calculator.g0                          # [m/s^2] Gravitational Constant
rho_0 = ISA_calculator.rho0                     # [kg/m^3] Air Density at Sea level
T0 = ISA_calculator.tempn[0]                    # [K] Temperature at Sea level

# Inputs from the collected data
Tm = Data_reader.flightdata['Dadc1_tat']                # [K] Measured Temperature
Vc = Data_reader.flightdata['Dadc1_tas']                # [m/s] Measured Airspeed                                               
hp = Data_reader.flightdata['Dadc1_alt']                # [N/m^2] Pressure altitude
delta_meas = Data_reader.flightdata['elevator_dte']     # [deg] Measured Elevator Deflection
AOA = Data_reader.flightdata['vane_AOA']                # [deg] Measured Angle of attack
Fe = Data_reader.flightdata['column_fe']                # [N] Measured elevator control force


Ve_bar_lst = np.zeros(len(Data_reader.flightdata['Dadc1_tas']))
red_el_def_lst = np.zeros(len(Data_reader.flightdata['elevator_dte']))
red_Fe_lst = np.zeros(len(Data_reader.flightdata['column_fe']))


def red_airspeed(hp, Vc, Tm):
    # alt = altitude
    # hp = measured pressure altitude 
    # Vc = calibrated airspeed on airplane airspeed indicator
    # Tm = measured temperature
    # Ve = equivalent airspeed
    # Ve_bar = reduced equivalent airspeed
    
    rho = ISA_calculator.ISAcalc(hp)[2]             # [kg/m^3]
                                       
    #calculation for the reduction of the measured airspeed
    p = p0 * (1 + (lam * hp)/T0) ** (g0/(lam * R))   
    M = m.sqrt(2/(gamma - 1) * ((1 + p0/p * ((1 + rho_0 * Vc ** 2 * (gamma - 1)/(2 * gamma * p0)) ** (gamma/(gamma -1)) -1)) ** ((gamma-1)/gamma) - 1))   
    T = (Tm + 273.15)/(1 + (gamma - 1)/2 * M ** 2)
    
    a = m.sqrt(gamma * R * T)   
    
    Vt = M * a
    Ve = Vt * m.sqrt(rho/rho_0)
    
    #calculation for the reduction of the equivalent airspeed
    m_tot = cg(100/2.2, payload_list)[1]
    W = m_tot * g0
    Ve_bar = Ve * m.sqrt(Ws/W)
    
    return Ve, Ve_bar


def red_thrust(delta_meas):
    #reduction of the non-standard engine thrust
    # delta_meas = measured elevator deflection
    # red_el_def = reduced elevator deflection
    
    Cmd = -1.1642       # [-] Elevator deflection moment coefficient
    Cmtc = -0.0064      # [-] Thrust moment arm
    Tcs = 0             # [-] Standard thrust coefficient
    Tc = 0              # [-] Thrust coefficient

    red_el_def = delta_meas - 1/Cmd * Cmtc * (Tcs -Tc)
    
    return red_el_def


def red_force(Fe):
    #calculations for the reduced elevator control force
    # Fe = measured elevator control force
    # red_Fe = reduced elevator control force
    
    #Aircraft weight calculations
    m_tot = cg(100/2.2, payload_list)[1]
    W = m_tot * g0
    
    red_Fe = Fe * Ws/W
    
    return red_Fe


def plotter(Ve_bar, el_def, F_e, alpha, CL, CD):
    #plotter for the different graphs
    
    fig = plt.figure()
    
    el_trim_curve = fig.add_subplot(221)
    el_trim_curve.set_title('Elevator Trim Curve')
    el_trim_curve.set_xlabel('Equivalent Airspeed')
    el_trim_curve.set_ylabel('Reduced Elevator Deflection')
    el_trim_curve.plot(Ve_bar, el_def)   
    
    el_force_curve = fig.add_subplot(222)
    el_force_curve.set_title('Reduced Elevator Control Force Curve')
    el_force_curve.set_xlabel('Equivalent Airspeed')
    el_force_curve.set_ylabel('Reduced Elevator Control Force')
    el_force_curve.plot(Ve_bar, F_e)
    
    Cl_alpha_curve = fig.add_subplot(223)
    Cl_alpha_curve.set_title('Cl vs Alpha Curve')
    Cl_alpha_curve.set_xlabel('Angle of Attack')
    Cl_alpha_curve.set_ylabel('Lift Coefficient')
    Cl_alpha_curve.plot(alpha, CL)
    
    CL_CD_curve = fig.add_subplot(224)
    CL_CD_curve.set_title('Cl vs Cd Curve')
    CL_CD_curve.set_xlabel('Drag Coefficient')
    CL_CD_curve.set_ylabel('Lift Coefficient')
    CL_CD_curve.plot(CD, CL)
    
    plt.show()


for i in range(0, len(Data_reader.flightdata['Dadc1_tas'])):
    Ve_bar_lst[i] = red_airspeed(hp[i], Vc[i], Tm[i])[1]
    red_el_def_lst[i] = red_thrust(delta_meas[i])
    red_Fe_lst[i] = red_force(Fe[i])



