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

alt = ft_to_m(flightdata['Dadc1_alt'])




#------------SYMM states-----------
for i in range(len(index)):
    index_begin = index[i][0]
    index_end = index[i][1]
    #mass calculation
    fuel_cons = pnd_to_kil(flightdata['rh_engine_FU'][index_begin]) + pnd_to_kil(flightdata['lh_engine_FU'][index_begin])
    mass_new = m-fuel_cons
    W_new = mass_new*g
    #---Startvals----
    start_alt = alt[index_begin]
    density_st, temper_st, press_st = ISAcalc(start_alt)
    vel_st = velocity[index_begin]
    alpha_st = alpha[index_begin]
    theta_st = theta[index_begin]
    q_st = q_val[index_begin]
    phi_st = phi[index_begin]
    p_st = p[index_begin]
    r_st = r[index_begin]


    muc = m/(density_st*S*c)
    mub = m/(density_st*S*b)
    CL_new = 2 * W_new / (density_st * vel_st ** 2 * S)
    CD_new = CD0 + (CLa * alpha_st ** 2) / (pi * A * e)
    CX0_new = W_new * sin(theta_st) / (0.5 * density_st * vel_st ** 2 * S)
    CZ0_new = -W_new*cos(theta_st)/(0.5*density_st*vel_st**2 * S)

    #Matrices for symmetric and asymmetric flight

    symmatc1 = symmetric_matrix_c1(muc,c, vel_st,CZadot,Cmadot,KY2)
    symmatc2 = symmetric_matrix_c2(CXu,CXa,CZ0_new,CZu,CZa,CX0,CZq,muc, Cmu,Cma, Cmq)
    symmatc3 = symmetric_matrix_c3(CXde,cxdt,CZde,czdt, Cmde,cmdt)
    asymmatc1 = asymmetric_matrix_c1(CYbdot,mub,b,vel_st,KX2,KXZ,Cnbdot,KZ2)
    asymmatc2 = asymmetric_matrix_c2(CYb,CL_new,CYp,CYr,mub, Clb, Clp, Clr,Cnb,Cnp,Cnr)
    asymmatc3 = asymmetric_matrix_c3(CYda,CYdr,Clda,Cldr,Cnda,Cndr)


    delta_e = deg_to_rad(flightdata['delta_e'])
    delta_e = delta_e[index_begin:index_end]
    trim = np.zeros(index_end-index_begin)
    ubar_symm = [delta_e,trim]

    delta_a = flightdata['delta_a'][index_begin:index_end]
    delta_r = flightdata['delta_r'][index_begin:index_end]

    ubar_asymm = [delta_a,delta_r]

    time_new = time[index_begin:index_end]

    abcd_sym = abcd_solver(symmatc1,symmatc2,symmatc3)
    abcd_asymm = abcd_solver(asymmatc1,asymmatc2,asymmatc3)
   # print(time_new)
    #print(len(time_new))
    ybar_symm = ybar(abcd_sym,ubar_symm,time_new)
    ybar_asymm = ybar(abcd_asymm,ubar_asymm,time_new)
    print('-------------------SYMMETRIC-----------------')
    print(ybar_symm)
    print('-------------------ASYMMETRIC-----------------')
    print(ybar_asymm)
    plt.plot(time_new,(vel_st+ybar_symm[:,0]*vel_st))
    plt.plot(time_new, velocity)
    plt.show()