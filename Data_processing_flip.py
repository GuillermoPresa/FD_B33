#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 15:04:57 2020

@author: frederik
"""

###### Data Processing tasks ########

import math as m
import cg_pos
import ISA_calculator

#### task 1: First Stationary Measurement Series

#inputs from the collected data
Tm = 0      # [K]
Vc = 0      # [m/s]
hp = 0      # [-]
alt = 0     # [m]


def red_airspeed(hp, Vc, Tm, alt):
    # alt = altitude
    # hp = measured pressure altitude 
    # Vc = calibrated airspeed on airplane airspeed indicator
    # Tm = measured temperature
    
    #input variables
    p0 = ISA_calculator.pressure0                   # [N/m^2]
    g0 = ISA_calculator.g0                          # [m/s^2]
    rho_0 = ISA_calculator.rho0                     # [kg/m^3]
    gamma = 1.4                                     # [-]
    R = ISA_calculator.R                            # [J/(mol*K)]
    lam = -0.0065                                   # [-]
    T0 = ISA_calculator.tempn[0]                    # [K]
    Actual_fuel_mass = cg_pos.Actual_fuel_mass      # [kg]
    rho = ISA_calculator.ISAcalc(alt)[2]            # [kg/m^3]
    Ws = 60500                                      # [N] from assingment
    
    #calculation for the reduction of the measured airspeed
    p = p0 * (1 + (lam * hp)/T0) ** (g0/(lam * R))   
    M = m.sqrt(2/(gamma - 1) * ((1 + p0/p * ((1 + rho_0 * Vc ** 2 * (gamma - 1)/(2 * gamma * p0)) ** (gamma/(gamma -1)) -1)) ** ((gamma-1)/gamma) - 1))   
    T = Tm/(1 + (gamma - 1)/2 * M ** 2)
    a = m.sqrt(gamma * R * T)
    Vt = M * a
    Ve = Vt * m.sqrt(rho/rho_0)
    
    #calculation for the reduction of the equivalent airspeed
    m_tot = cg_pos.cg(Actual_fuel_mass)[1]
    W = m_tot * g0
    Ve_bar = Ve * m.sqrt(Ws/W)
    
    return Ve, Ve_bar


def red_engine_thrust():    
    #reduction of the non-standard engine thrust
    
    
