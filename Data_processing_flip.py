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
import Main

#### task 1: First Stationary Measurement Series

# Aircraft inputs
S = 30.00                                               # [m^2] Wing Surface Area
c_bar = 2.0569                                          # [m] Average Chord Length

# Envirnonment inputs
gamma = 1.4                                             # [-]
lam = -0.0065                                           # [-]
Ws = 60500                                              # [N] Reduced aircraft weight from assingment
R = 287.0                                               # [J/(mol*K)] Specific gas Constant of Air

# Sea level conditions
p0 = ISA_calculator.pressure0                           # [N/m^2] Pressure at Sea level
g0 = ISA_calculator.g0                                  # [m/s^2] Gravitational Constant
rho_0 = ISA_calculator.rho0                             # [kg/m^3] Air Density at Sea level
T0 = ISA_calculator.tempn[0]                            # [K] Temperature at Sea level

# Inputs from the collected data
Tm = Data_reader.flightdata['Dadc1_tat']                # [K] Measured Temperature
Vc = Data_reader.flightdata['Dadc1_tas']                # [m/s] Measured Airspeed                                               
hp = Data_reader.flightdata['Dadc1_alt']                # [N/m^2] Pressure altitude
el_def_meas = Data_reader.flightdata['elevator_dte']     # [deg] Measured Elevator Deflection
AOA = Data_reader.flightdata['vane_AOA']                # [deg] Measured Angle of attack
Fe = Data_reader.flightdata['column_fe']                # [N] Measured elevator control force
time = Data_reader.flightdata['time']                   # [s] Time

#to be calculated
AFM = np.zeros(len(Data_reader.flightdata['Dadc1_tas'])) #actual fuel mass till needs to be calculated 
T = np.zeros(len(Data_reader.flightdata['Dadc1_tas'])) #thrust still needs to be calculated from thrust.exe


Ve_bar_lst = np.zeros(len(Data_reader.flightdata['Dadc1_tas']))
red_el_def_lst = np.zeros(len(Data_reader.flightdata['elevator_dte']))
red_Fe_lst = np.zeros(len(Data_reader.flightdata['column_fe']))


def red_airspeed(hp, Vc, Tm, AFM):
    # alt = altitude
    # hp = measured pressure altitude 
    # Vc = calibrated airspeed on airplane airspeed indicator
    # Tm = measured temperature
    # Ve = equivalent airspeed
    # Ve_bar = reduced equivalent airspeed
    # AFM = actual fuel mass
    
    rho = ISA_calculator.ISAcalc(hp)[2]             # [kg/m^3]
                                       
    #calculation for the reduction of the measured airspeed
    p = p0 * (1 + (lam * hp)/T0) ** (g0/(lam * R))   
    M = m.sqrt(2/(gamma - 1) * ((1 + p0/p * ((1 + rho_0 * Vc ** 2 * (gamma - 1)/(2 * gamma * p0)) ** (gamma/(gamma -1)) -1)) ** ((gamma-1)/gamma) - 1))   
    T = (Tm + 273.15)/(1 + (gamma - 1)/2 * M ** 2)
    
    a = m.sqrt(gamma * R * T)   
    
    Vt = M * a
    Ve = Vt * m.sqrt(rho/rho_0)
    
    #calculation for the reduction of the equivalent airspeed
    m_tot = cg(AFM, payload_list)[1]
    W = m_tot * g0
    Ve_bar = Ve * m.sqrt(Ws/W)
    
    return Ve, Ve_bar, rho, W, T


def red_thrust(el_def_meas):
    #reduction of the non-standard engine thrust
    # delta_meas = measured elevator deflection
    # red_el_def = reduced elevator deflection
    
    Cmd = -1.1642       # [-] Elevator deflection moment coefficient
    Cmtc = -0.0064      # [-] Thrust moment arm
    Tcs = 0             # [-] Standard thrust coefficient
    Tc = 0              # [-] Thrust coefficient

    red_el_def = el_def_meas - 1/Cmd * Cmtc * (Tcs -Tc)
    
    return red_el_def


def red_force(Fe, hp, Vc, Tm, AFM):
    #calculations for the reduced elevator control force
    # Fe = measured elevator control force
    # red_Fe = reduced elevator control force
    
    red_Fe = Fe * Ws/(red_airspeed(hp, Vc, Tm, AFM)[3])
    
    return red_Fe


def C_md(time1, time2):
    #calculations for the determination of the elevator effectiveness coefficient C_md
    # time1 = time at beginning of weight shift (Xcg1)
    # time2 = time at end of weight shift (Xcg2)
    # el_def1 = elevator deflection at beginning of weight shift
    # el_def2 = elevator deflection at end of weight shift
    # C_md = Elevator effectiveness coefficient
    
    abs_difference1 = lambda list_value : abs(list_value - time1)
    closest1 = min(time, key=abs_difference1)
    index1 = time.index(closest1)
    
    abs_difference2 = lambda list_value : abs(list_value - time2)
    closest2 = min(time, key=abs_difference2) 
    index2 = time.index(closest2)
    
    el_def1 = el_def_meas[index1]
    el_def2 = el_def_meas[index2]
    
    C_md = - 1/(el_def2 - el_def1) * Main.CN * (Main.Xcg2 - Main.Xcg1)/c_bar
    
    return C_md


def plotter_stat_meas1(alpha, CL, CD, CL2):
    #plotter for the first static measurement graphs
    
    fig = plt.figure()
    
    
    Cl_alpha_curve = fig.add_subplot(221)
    Cl_alpha_curve.set_title('Cl vs Alpha Curve')
    Cl_alpha_curve.set_xlabel('Angle of Attack')
    Cl_alpha_curve.set_ylabel('Lift Coefficient')
    Cl_alpha_curve.plot(alpha, CL)
    
    Cd_alpha_curve = fig.add_subplot(222)
    Cd_alpha_curve.set_title('Cd vs Alpha Curve')
    Cd_alpha_curve.set_xlabel('Angle of Attack')
    Cd_alpha_curve.set_ylabel('Drag Coefficient')
    Cd_alpha_curve.plot(alpha, CD)
    
    CL_CD_curve = fig.add_subplot(223)
    CL_CD_curve.set_title('Cl vs Cd Curve')
    CL_CD_curve.set_xlabel('Drag Coefficient')
    CL_CD_curve.set_ylabel('Lift Coefficient')
    CL_CD_curve.plot(CD, CL)
    
    CD_CL2_curve = fig.add_subplot(223)
    CD_CL2_curve.set_title('Cd vs Cl^2 Curve')
    CD_CL2_curve.set_xlabel('Drag Coefficient')
    CD_CL2_curve.set_ylabel('Squared Lift Coefficient')
    CD_CL2_curve.plot(CD, CL2)
    
    plt.show()
    
    
def plotter_stat_meas2(Ve_bar, el_def, F_e):
    #plotter for the second static measurement graphs
    
    fig = plt.figure()
    
    el_trim_curve = fig.add_subplot(211)
    el_trim_curve.set_title('Elevator Trim Curve')
    el_trim_curve.set_xlabel('Equivalent Airspeed')
    el_trim_curve.set_ylabel('Reduced Elevator Deflection')
    el_trim_curve.plot(Ve_bar, el_def)   
    
    el_force_curve = fig.add_subplot(212)
    el_force_curve.set_title('Reduced Elevator Control Force Curve')
    el_force_curve.set_xlabel('Equivalent Airspeed')
    el_force_curve.set_ylabel('Reduced Elevator Control Force')
    el_force_curve.plot(Ve_bar, F_e)
    
    plt.show()


############## Calculations with main.py functions

# substitute reduced values in static measuerements
Static_Measurements_1 = Main.Static_Measurements_1
Static_Measurements_2 = Main.Static_Measurements_2
Static_Measurements_3 = Main.Static_Measurements_3

Static_Measurements_1.DataLineList = Main.Static_Measurements_1.DataLineList
Static_Measurements_2.DataLineList = Main.Static_Measurements_2.DataLineList
Static_Measurements_3.DataLineList = Main.Static_Measurements_3.DataLineList


for j in range(0, len(Main.Static_Measurements_1.DataLineList)):
    Static_Measurements_1.DataLineList[j][1] = red_airspeed(Main.Static_Measurements_1.DataLineList[j][0], Main.Static_Measurements_1.DataLineList[j][1], Main.Static_Measurements_1.DataLineList[j][6], (Main.TotalFuelMass - Main.Static_Measurements_1.DataLineList[j][5]))[1]
    Static_Measurements_1.DataLineList[j][6] = red_airspeed(Main.Static_Measurements_1.DataLineList[j][0], Main.Static_Measurements_1.DataLineList[j][1], Main.Static_Measurements_1.DataLineList[j][6], (Main.TotalFuelMass - Main.Static_Measurements_1.DataLineList[j][5]))[4]


for j in range(0, len(Main.Static_Measurements_2.DataLineList)):
    Static_Measurements_2.DataLineList[j][1] = red_airspeed(Main.Static_Measurements_2.DataLineList[j][0], Main.Static_Measurements_2.DataLineList[j][1], Main.Static_Measurements_2.DataLineList[j][10], (Main.TotalFuelMass - Main.Static_Measurements_2.DataLineList[j][8]))[1]
    Static_Measurements_2.DataLineList[j][10] = red_airspeed(Main.Static_Measurements_2.DataLineList[j][0], Main.Static_Measurements_2.DataLineList[j][10], Main.Static_Measurements_2.DataLineList[j][10], (Main.TotalFuelMass - Main.Static_Measurements_2.DataLineList[j][8]))[4]


for j in range(0, len(Main.Static_Measurements_3.DataLineList)):
    Static_Measurements_3.DataLineList[j][1] = red_airspeed(Main.Static_Measurements_3.DataLineList[j][0], Main.Static_Measurements_3.DataLineList[j][1], Main.Static_Measurements_3.DataLineList[j][10], (Main.TotalFuelMass - Main.Static_Measurements_3.DataLineList[j][8]))[1]
    Static_Measurements_3.DataLineList[j][10] = red_airspeed(Main.Static_Measurements_3.DataLineList[j][0], Main.Static_Measurements_3.DataLineList[j][1], Main.Static_Measurements_3.DataLineList[j][10], (Main.TotalFuelMass - Main.Static_Measurements_3.DataLineList[j][8]))[4]


# plotting the values for stat meas 1
alpha = Main.Coeficients1(Static_Measurements_1)[0]
CL = Main.Coeficients1(Static_Measurements_1)[1]
CD = Main.Coeficients1(Static_Measurements_1)[2]
CL2 = np.square(CL)

plotter_stat_meas1(alpha, CL, CD, CL2)
print(CL)
print(CD)

# plotting the values for stat meas 2


# for i in range(0, len(Data_reader.flightdata['Dadc1_tas'])):
#     Ve_bar_lst[i] = red_airspeed(hp[i], Vc[i], Tm[i], AFM[i])[1]
#     red_el_def_lst[i] = red_thrust(el_def_meas[i])
#     red_Fe_lst[i] = red_force(Fe[i], hp[i], Vc[i], Tm[i], AFM[i])



