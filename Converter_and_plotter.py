from Cit_par import *
from Data_reader import *
import matplotlib.pyplot as plt
from math import *
from state_space import *
phu_index = 0
spir_index = 0
short_index = 0
dutch_index = 0
aproll_index = 0
t = 0


def psi_to_pascal(psi_vals):
    return psi_vals*6894.76
def knots_to_mps(knots_vals):
    return knots_vals*0.514444
def deg_to_rad(deg_vals):
    return deg_vals*pi/180
W_new = m*g
rho_new = rho0 * pow(((1+(lambda1 * hp0 / Temp0))), (-((g / (lambda1*R)) + 1)))
muc = m/(rho_new*S*b)
CL_new = 2*W/(rho_new*V0**2 * S)
CD_new = CD0 + (CLa*alpha0**2)/(pi*A*e)
CX0_new = W*sin(th0)/(0.5*rho_new*V0**2 * S)