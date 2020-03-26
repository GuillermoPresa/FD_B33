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
index = [(32300,33300),(39000,39150),(34000,34500),(35130,35300),(37110,38000)]
#index = [(31900,33200),(38800,39800),(34210,34300),(35130,35300),(37000,37200)]
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
yaw_beta = [0]
phi = deg_to_rad(flightdata['Ahrs1_Roll'])
p = deg_to_rad(flightdata['Ahrs1_bRollRate'])
r = deg_to_rad(flightdata['Ahrs1_bYawRate'])
for i in range(len(r)-1):
    yaw_beta.append(yaw_beta[i] + r[i])
alt = ft_to_m(flightdata['Dadc1_alt'])
yaw_beta = np.asarray(yaw_beta)
def error(real,sim):
    return sqrt(sum(((real-sim)/real)**2))/len(real)
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
    yaw_beta_st = yaw_beta[index_begin]


    muc = m/(density_st*S*c)
    mub = m/(density_st*S*b)
    CL_new = 2 * W_new / (density_st * (vel_st ** 2) * S)
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

    #matasymeig = np.linalg.eigvals(-np.matmul(np.linalg.inv(symmatc1),symmatc2))
    #mataasymeig = np.linalg.eigvals(-np.matmul(np.linalg.inv(asymmatc1),asymmatc2))
    delta_e = deg_to_rad(flightdata['delta_e'])
    delta_e = delta_e[index_begin:index_end]
    trim = np.zeros(index_end-index_begin)
    ubar_symm = delta_e

    delta_a = flightdata['delta_a'][index_begin:index_end]
    delta_r = flightdata['delta_r'][index_begin:index_end]

    ubar_asymm = [delta_a,delta_r]

    time_new = time[index_begin:index_end]
   # tns = time_new
   # tna = time_new
    print(manuever)
    abcd_symm = symabcd_solver(symmatc1,symmatc2,symmatc3)
    abcd_asymm = asymabcd_solver(asymmatc1,asymmatc2,asymmatc3)
   # print(time_new)
   # time_new = np.linspace(time(index_begin),time(index_end),index_end-index_begin)
    ybar_symm,tns,x1 = ybar(abcd_symm,ubar_symm,time_new,[0,0,theta_st,q_st]) #vel_st,alpha_st,theta_st,q_st]
    ybar_asymm,tna,x2 = ybar(abcd_asymm,ubar_asymm,time_new,[0,phi_st,p_st,r_st])#[yaw_beta_st,phi_st,p_st,q_st]








    #print(ybar_asymm)

    if manuever_type[i] == 1:
        print('-------------------SYMMETRIC-----------------')
       #  plt.figure(figsize =(16,12))
       #  plt.subplot(2,1,1)
       #  plt.plot(tns, ybar_symm[0]*vel_st + vel_st, label = 'Velocity [m/s]')
       #  #plt.plot(tns,ybar_symm[1]+alpha_st, label = 'Alpha [Rad]')
       #  #plt.plot(tns,ybar_symm[2]+theta_st, label = 'Theta [Rad]')
       # # plt.plot(tns,ybar_symm[3]*vel_st/c, label = 'q [Rad/s]')
       #  plt.xlabel('Time Steps [s]')
       #  plt.ylabel(' Velocity m/s')
       #  plt.subplot(2,1,2)
       #  plt.plot(tns,deg_to_rad(flightdata['delta_e'][index_begin:index_end]), label = 'delta_e')
       #  plt.xlabel('Time Steps [s]')
       #  plt.ylabel('delta e [Rad]')
       #  plt.show()
        plt.figure(figsize=(14, 10))
        plt.subplot(2,2,1)
        plt.plot(tns,(ybar_symm[0]*vel_st + vel_st), label = 'Simulated Response')
        plt.plot(tns, velocity[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time Steps')
        plt.ylabel('Velocity [m/s]')
        plt.subplot(2,2,2)
        plt.plot(tns,ybar_symm[1]+alpha_st, label = 'Simulated Response')
        plt.plot(tns,alpha[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time Steps')
        plt.ylabel('Alpha [Rad]')
        plt.subplot(2,2,3)
        plt.plot(tns, ybar_symm[2]+theta_st, label = 'Simulated Response')
        plt.plot(tns,theta[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time Steps')
        plt.ylabel('Theta [Rad]')
        plt.subplot(2,2,4)
        plt.plot(tns, ybar_symm[3]*vel_st/c, label = 'Simulated Response')
        plt.plot(tns, q_val[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time Steps')
        plt.ylabel('Pitch Rate [Rad/s]')
        plt.suptitle(manuever)
        plt.show()
        av_error = error(velocity[index_begin:index_end],ybar_symm[0]*vel_st + vel_st),error(alpha[index_begin:index_end],ybar_symm[1]+alpha_st), error(theta[index_begin:index_end],ybar_symm[2]+theta_st),error(q_val[index_begin:index_end],ybar_symm[3]*vel_st/c)
        print(av_error)
    else:
    #Asymm Plots
        print('-------------------ASYMMETRIC-----------------')
       # print(mataasymeig)
        plt.figure(figsize=(16,12))
        plt.subplot(2,2,1)
        plt.plot(tna,ybar_asymm[0]+yaw_beta_st, label = 'Simulated Response')
        plt.plot(tna, yaw_beta[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time Steps')
        plt.ylabel('Beta [Rad]')
        plt.subplot(2,2,2)
        plt.plot(tna,ybar_asymm[1], label = 'Simulated Response')
        plt.plot(tna,phi[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time Steps')
        plt.ylabel('Phi [Rad]')
        plt.subplot(2,2,3)
        plt.plot(tna, -ybar_asymm[2], label = 'Simulated Response')
        plt.plot(tna,p[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time Steps')
        plt.ylabel('p Values [Rad/s]')
        plt.subplot(2,2,4)
        plt.plot(tna, ybar_asymm[3]+r_st, label = 'Simulated Response')
        plt.plot(tna, r[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time Steps')
        plt.ylabel('r Values [Rad/s]')
        plt.suptitle(manuever)
        plt.show()
        av_error = error(yaw_beta[index_begin:index_end],ybar_asymm[0]+yaw_beta_st), error(phi[index_begin:index_end],ybar_asymm[1]), error(p[index_begin:index_end],-ybar_asymm[2]), error(r[index_begin:index_end],ybar_asymm[3]+r_st)
        print(av_error)
      #   plt.figure(figsize = (16,12))
      #   plt.subplot(2,1,1)
      # #  plt.plot(tna,ybar_asymm[0]+yaw_beta_st, label = 'Beta [Rad]')
      #   plt.plot(tna,ybar_asymm[1],label = 'Phi [Rad]')
      #  # plt.plot(tna,-ybar_asymm[2],label = 'p [Rad/s]')
      #   #plt.plot(tna, ybar_asymm[3]+r_st, label = 'r [Rad/s]')
      #   #plt.legend(loc = 'best')
      #   plt.xlabel('Time Steps [s]')
      #   plt.ylabel('Phi [Rad]')
      #   plt.subplot(2,1,2)
      #   plt.plot(tna,deg_to_rad(flightdata['delta_r'][index_begin:index_end]),label = 'delta r')
      #   plt.plot(tna,deg_to_rad(flightdata['delta_a'][index_begin:index_end]), label = 'delta a')
      #   plt.legend(loc = 'best')
      #   plt.xlabel('Time Steps[s]')
      #   plt.ylabel('delta a, delta r [Rad]')
      #   plt.show()

