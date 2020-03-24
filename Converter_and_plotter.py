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
index = [(32470,33970),(34210,35410),(35130,35330),(36000,36150),(37110,37260)]
name = ['Phugoid', 'Spiral', 'Short Period', 'Dutch Roll', 'Aperiodic']
#Manuever type 1 means symmetric and 0 means asymmetric
manuever_type = [1,0,1,0,0]
time = flightdata['time'] - 9
#----SYMM stuff---
velocity = knots_to_mps(flightdata['Dadc1_tas'])
alpha = deg_to_rad(flightdata['vane_AOA'])
theta = deg_to_rad(flightdata['Ahrs1_Pitch'])
q_val = deg_to_rad(flightdata['Ahrs1_bPitchRate'])
#---ASYMM stuff---
yaw_beta = deg_to_rad(flightdata['Ahrs1_bYawRate'])
phi = deg_to_rad(flightdata['Ahrs1_Roll'])
p = deg_to_rad(flightdata['Ahrs1_bRollRate'])
r = deg_to_rad(flightdata['Ahrs1_bYawRate'])

alt = ft_to_m(flightdata['Dadc1_alt'])




#------------SYMM states-----------
for i in range(len(index)):

    #These define the start and the end points of an index
    index_begin = index[i][0]
    index_end = index[i][1]

    manuever = name[i]

    #mass calculation
    fuel_cons = pnd_to_kil(flightdata['rh_engine_FU'][index_begin]) + pnd_to_kil(flightdata['lh_engine_FU'][index_begin])
    mass_new = m-fuel_cons
    W_new = mass_new*g
    #---Startvals----
    #Values at the beginning of each index to be used in the calculations below
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

    matasymeig = np.linalg.eigvals(-np.matmul(np.linalg.inv(symmatc1),symmatc2))
    mataasymeig = np.linalg.eigvals(np.matmul(np.linalg.inv(asymmatc1),asymmatc2))
    delta_e = deg_to_rad(flightdata['delta_e'])
    delta_e = delta_e[index_begin:index_end]
    trim = np.zeros(index_end-index_begin)
    ubar_symm = [delta_e,trim]

    delta_a = flightdata['delta_a'][index_begin:index_end]
    delta_r = flightdata['delta_r'][index_begin:index_end]

    ubar_asymm = [delta_a,delta_r]

    time_new = time[index_begin:index_end]

    abcd_sym = symabcd_solver(symmatc1,symmatc2,symmatc3)
    abcd_asymm = asymabcd_solver(asymmatc1,asymmatc2,asymmatc3)
   # print(time_new)

    ybar_symm = ybar(abcd_sym,ubar_symm,time_new,np.array([0,0,theta_st,q_st]))
    ybar_asymm = ybar(abcd_asymm,ubar_asymm,time_new,np.array([[0],[phi_st],[p_st],[r_st]]))
    print(manuever)







    #print(ybar_asymm)

    if manuever_type[i] == 1:
        print('-------------------SYMMETRIC-----------------')
        print(matasymeig)
        plt.figure(figsize=(14, 10))
        plt.subplot(2,2,1)
        plt.plot(time_new,(vel_st+ybar_symm[0]*vel_st), label = 'Simulated Response')
        plt.plot(time_new, velocity[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time Steps')
        plt.ylabel('Velocity [m/s]')
        plt.subplot(2,2,2)
        plt.plot(time_new,ybar_symm[1]+alpha_st, label = 'Simulated Response')
        plt.plot(time_new,alpha[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time Steps')
        plt.ylabel('Alpha [Rad]')
        plt.subplot(2,2,3)
        plt.plot(time_new, ybar_symm[2]+theta_st, label = 'Simulated Response')
        plt.plot(time_new,theta[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time Steps')
        plt.ylabel('Theta [Rad]')
        plt.subplot(2,2,4)
        plt.plot(time_new, ybar_symm[3]*(vel_st/c), label = 'Simulated Response')
        plt.plot(time_new, q_val[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time Steps')
        plt.ylabel('Pitch Rate [Rad/s]')
        plt.suptitle(manuever)
        plt.show()
    else:
    #Asymm Plots
        print('-------------------ASYMMETRIC-----------------')
        print(mataasymeig)
        plt.figure(figsize=(16,12))
        plt.subplot(2,2,1)
        plt.plot(time_new,(ybar_asymm[0]), label = 'Simulated Response')
        plt.plot(time_new, yaw_beta[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time Steps')
        plt.ylabel('Beta [Rad]')
        plt.subplot(2,2,2)
        plt.plot(time_new,phi_st+ybar_asymm[1], label = 'Simulated Response')
        plt.plot(time_new,phi[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time Steps')
        plt.ylabel('Phi [Rad]')
        plt.subplot(2,2,3)
        plt.plot(time_new, (ybar_asymm[2]*2*vel_st/b)-p_st, label = 'Simulated Response')
        plt.plot(time_new,p[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time Steps')
        plt.ylabel('p Values [Rad/s]')
        plt.subplot(2,2,4)
        plt.plot(time_new, (ybar_asymm[3]*2*vel_st/b)+r_st, label = 'Simulated Response')
        plt.plot(time_new, r[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time Steps')
        plt.ylabel('r Values [Rad/s]')
        plt.suptitle(manuever)
        plt.show()

