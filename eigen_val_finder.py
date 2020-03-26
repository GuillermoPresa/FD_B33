from Cit_par import *
from Data_reader import *
import matplotlib.pyplot as plt
from math import *
from ISA_calculator import ISAcalc
from state_space import *
#Equations to convert units
def psi_to_pascal(psi_vals):
    return psi_vals*6894.76
def knots_to_mps(knots_vals):
    return knots_vals*0.514444
def deg_to_rad(deg_vals):
    return deg_vals*pi/180
def ft_to_m(h_vals):
    return h_vals*0.3048
def pnd_to_kil(w_vals):
    return 0.453592*w_vals
index = [(32300,33300),(39000,39150),(34000,34500),(35130,35300),(37110,38000)]
name = ['Phugoid', 'Spiral', 'Short Period', 'Dutch Roll', 'Aperiodic']
velocity_phug = knots_to_mps(flightdata['Dadc1_tas'][32300:33300])
theta_short = deg_to_rad(flightdata['Ahrs1_Pitch'][39000:39150])
time =  flightdata['time'] - 9
f1,f2,f3,f4 = np.polyfit(time[32300:33300],velocity_phug,3)
g1,g2,g3,g4 = np.polyfit(time[39000:39150],theta_short,3)
def func(val,a,b,c,d):
    return a*val**3 + b*val**2 + c*val + d
plt.plot(time[32300:33300],velocity_phug)
plt.plot(np.linspace(0,1000,2000),func(np.linspace(0,1000,2000),f1,f2,f3,f4))
plt.show()