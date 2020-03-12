from Cit_par import *
from Data_reader import *
import matplotlib.pyplot as plt
from math import *
from state_space import *

def psi_to_pascal(psi_vals):
    return psi_vals*6894.76
def knots_to_mps(knots_vals):
    return knots_vals*0.514444
def deg_to_rad(deg_vals):
    return deg_vals*pi/180

