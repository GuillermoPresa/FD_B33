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
Tm = 0
Vc = 0
hp = 0
alt = 0
rho = ISA_calculator.ISAcalc(alt)[2]

#standard inputs
p0 = ISA_calculator.pressure0
g0 = ISA_calculator.g0
rho_0 = ISA_calculator.rho0
gamma = 1.4
R = ISA_calculator.R
lam = -0.0065
T0 = ISA_calculator.tempn[0]
Actual_fuel_mass = cg_pos.Actual_fuel_mass


def red_measured_airspeed(p0, lam, hp, T0, g0, R, gamma, Vc, Tm, rho_0, rho):
    #calculation for the reduction of the measured airspeed
    # p0 = pressure at sea level
    # lam = lambda
    # hp = measured pressure altitude 
    # T0 = temperature at sea level
    # g0 = 
    # R = specific gas constant of air
    # gamma
    # Vc = calibrated airspeed on airplane airspeed indicator
    # Tm = measured temperature
    # rho_0 = air density at sea level
    # rho = measured air density
    
    p = p0 * (1 + (lam * hp)/T0) ** (g0/(lam * R))   
    M = m.sqrt(2/(gamma - 1) * ((1 + p0/p * ((1 + rho_0 * Vc ** 2 * (gamma - 1)/(2 * gamma * p0)) ** (gamma/(gamma -1)) -1)) ** ((gamma-1)/gamma) - 1))   
    T = Tm/(1 + (gamma - 1)/2 * M ** 2)
    a = m.sqrt(gamma * R * T)
    Vt = M * a
    Ve = Vt * m.sqrt(rho/rho_0)
    
    return Ve


m_tot = cg_pos.cg(Actual_fuel_mass)[1]
W = m_tot * g0
# Ws = 

Ve_bar = red_measured_airspeed(p0, lam, hp, T0, g0, R, gamma, Vc, Tm, rho_0, rho) * m.sqrt(Ws/W)
    


#reduction of the non-standard engine thrust

def red_engine_thrust():    
    
    
    
