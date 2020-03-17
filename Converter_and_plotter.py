from Cit_par import *
from Data_reader import *
import matplotlib.pyplot as plt
from math import *
from state_space import *
#Equations to convert units
def psi_to_pascal(psi_vals):
    return psi_vals*6894.76
def knots_to_mps(knots_vals):
    return knots_vals*0.514444
def deg_to_rad(deg_vals):
    return deg_vals*pi/180
#phugoid, spiral, short period, dutch roll, aperiodic indexes
index = [(32470,34210),(34210,35130),(35130,36000),(36000,37110),(37110,39000)]
time = flightdata['time'] - 9
#----SYMM stuff---
velocity = knots_to_mps(flightdata['Dadc1_tas'])
alpha = deg_to_rad(flightdata['vane_AOA'])
theta = deg_to_rad(flightdata['Ahrs1_Pitch'])
q_val = deg_to_rad(flightdata['Ahrs1_bPitchRate'])
#---ASYMM stuff---
phi = deg_to_rad(flightdata['Ahrs1_Roll'])
p = deg_to_rad(flightdata['Ahrs1_bRollRate'])
r = deg_to_rad(flightdata['Ahrs1_bYawRate'])



#------------SYMM states-----------
for i in range(len(index)):
    index_begin = index[i][0]
    index_end = index[i][1]



W_new = m*g
rho_new = rho0 * pow(((1+(lambda1 * hp0 / Temp0))), (-((g / (lambda1*R)) + 1)))
muc = m/(rho_new*S*b)
CL_new = 2*W/(rho_new*V0**2 * S)
CD_new = CD0 + (CLa*alpha0**2)/(pi*A*e)
CX0_new = W*sin(th0)/(0.5*rho_new*V0**2 * S)


