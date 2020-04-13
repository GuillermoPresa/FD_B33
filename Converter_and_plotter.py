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
#These are the time indexes used to find the start and end points of each eigenmotion
#They may need to be adjusted as these are based on observations in the recorded data, as the
#data given within the excel sheet is faulty

if file == 'FTISxprt-20200306_flight1.mat':
    #new index finder as from the flight data
    t_begin_ph = 3195       #begin time phugoid in seconds 
    t_begin_sp = 3418       #begin time short period in seconds
    t_begin_dr = 3513       #begin time dutch roll in seconds 
    t_begin_ar = 3710       #begin time aperiodic roll in seconds 
    t_begin_spir = 3880     #begin time spiral in seconds 
else:
    #refernce data
    t_begin_ph = 53*60 + 57       #begin time phugoid in seconds 
    t_begin_sp = 60*60 + 35       #begin time short period in seconds
    t_begin_dr = 61*60 + 57       #begin time dutch roll in seconds
    t_begin_ar = 59*60 + 10       #begin time aperiodic roll in seconds 
    t_begin_spir = 65*60 + 20     #begin time spiral in seconds

#approximate length maneuvre [s]
len_ph = 200
len_sp = 10
len_dr = 15
len_ar = 15
len_spir = 150


#indexes
i_ph = int(np.where(flightdata['time']==t_begin_ph)[0])
i_sp = int(np.where(flightdata['time']==t_begin_sp)[0])
i_dr = int(np.where(flightdata['time']==t_begin_dr)[0])
i_ar = int(np.where(flightdata['time']==t_begin_ar)[0])
i_spir = int(np.where(flightdata['time']==t_begin_spir)[0])

index = [(i_ph, i_ph+len_ph*10), (i_spir, i_spir+len_spir*10), (i_sp, i_sp+len_sp*10), (i_dr, i_dr+len_dr*10), (i_ar, i_ar+len_ar*10)]

#index = [(32300,33300),(39000,39150),(34000,34500),(35130,35300),(37110,38000)]
#index = [(31900,33200),(38800,39800),(34210,34300),(35130,35300),(37000,37200)]


name = ['Phugoid', 'Spiral', 'Short Period', 'Dutch Roll', 'Aperiodic Roll']
#Manuever type 1 means symmetric and 0 means asymmetric
manuever_type = [1,0,1,0,0]
#This line below resets the time to start at 0 because the time started at 9
time = flightdata['time']


#This is the data required to plot and compare the symmetric values
velocity = knots_to_mps(flightdata['Dadc1_tas'])
alpha = deg_to_rad(flightdata['vane_AOA'])
theta = deg_to_rad(flightdata['Ahrs1_Pitch'])
q_val = deg_to_rad(flightdata['Ahrs1_bPitchRate'])


#This is the data required to plot and compare the asymmetric values. Since
#there were no yaw values provided, the yaw values had to be calculated using a loop and the yaw rate (r)
yaw_beta = [0]
phi = deg_to_rad(flightdata['Ahrs1_Roll'])
p = deg_to_rad(flightdata['Ahrs1_bRollRate'])
r = deg_to_rad(flightdata['Ahrs1_bYawRate'])
for i in range(len(r)-1):
    #yaw_beta.append(yaw_beta[i] + r[i])
    yaw_beta.append(yaw_beta[i])
alt = ft_to_m(flightdata['Dadc1_alt'])
yaw_beta = np.asarray(yaw_beta)

#This can be used to calculate the error in the simulated values versus the real values
def tot_error(real,sim):
    return sqrt(sum((real-sim)**2))

def RMS(real,sim):
    RMS = np.sqrt(np.mean((real-sim)**2))
    return RMS

tot_error_sym = []
tot_error_asym = []
    
#------------SYMM states-----------
for i in range(len(index)):

    #These define the start and the end points of an index
    index_begin = index[i][0]
    index_end = index[i][1]

    manuever = name[i]

    #mass calculation this was used to factor in the changing weight of the aircraft due to the mass of fuel changing
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

    #These are the recalculated values that changed with each manoeuver
    #In order to change the CLa, e and Cma values go into the cit_par.py file and change them there
    muc = m/(density_st*S*c)
    mub = m/(density_st*S*b)
    CL_new = 2 * W_new / (density_st * (vel_st ** 2) * S)
    CD_new = CD0 + (CLa * alpha_st ** 2) / (pi * A * e)
    CX0_new = W_new * sin(theta_st) / (0.5 * density_st * vel_st ** 2 * S)
    CZ0_new = -W_new*cos(theta_st)/(0.5*density_st*vel_st**2 * S)

    #Matrices for symmetric and asymmetric flight
    #These matrices come from the state_space.py import

    symmatc1 = symmetric_matrix_c1(muc,c, vel_st,CZadot,Cmadot,KY2)
    symmatc2 = symmetric_matrix_c2(CXu,CXa,CZ0_new,CZu,CZa,CX0,CZq,muc, Cmu,Cma, Cmq)
    # symmatc3 = symmetric_matrix_c3(CXde,cxdt,CZde,czdt, Cmde,cmdt)
    symmatc3 = symmetric_matrix_c3(CXde,CZde,Cmde)
    asymmatc1 = asymmetric_matrix_c1(CYbdot,mub,b,vel_st,KX2,KXZ,Cnbdot,KZ2)
    asymmatc2 = asymmetric_matrix_c2(CYb,CL_new,CYp,CYr,mub, Clb, Clp, Clr,Cnb,Cnp,Cnr)
    asymmatc3 = asymmetric_matrix_c3(CYda,CYdr,Clda,Cldr,Cnda,Cndr)

    #matasymeig = np.linalg.eigvals(-np.matmul(np.linalg.inv(symmatc1),symmatc2))
    #mataasymeig = np.linalg.eigvals(-np.matmul(np.linalg.inv(asymmatc1),asymmatc2))
    delta_e = deg_to_rad(flightdata['delta_e'])
    delta_e = delta_e[index_begin:index_end]
    trim = np.zeros(index_end-index_begin)
    ubar_symm = delta_e

    delta_a = deg_to_rad(flightdata['delta_a'][index_begin:index_end])
    delta_r = deg_to_rad(flightdata['delta_r'][index_begin:index_end])
    
    ubar_asymm = np.vstack((delta_a, delta_r))

    # time_new = time[index_begin:index_end]
    time_new = np.linspace(0, (index_end - index_begin)/10, (index_end - index_begin))
    
   # tns = time_new
   # tna = time_new
    print('-------------------',manuever,'-------------------')
    if manuever_type[i] == 1:
        abcd_symm = symabcd_solver(symmatc1,symmatc2,symmatc3)
        tns,ybar_symm,x1 = ybar(abcd_symm,ubar_symm,time_new) #,[0,0,theta_st,q_st]) #vel_st,alpha_st,theta_st,q_st]
    else:
        abcd_asymm = asymabcd_solver(asymmatc1,asymmatc2,asymmatc3)
        tna,ybar_asymm,x2 = ybar(abcd_asymm,ubar_asymm, time_new) #,[0,phi_st,p_st,r_st])#[yaw_beta_st,phi_st,p_st,q_st]
   # print(time_new)
   # time_new = np.linspace(time(index_begin),time(index_end),index_end-index_begin)
    
    #print(ybar_asymm)

    if manuever_type[i] == 1:
        print('-------------------SYMMETRIC-----------------')
        # plt.figure(figsize =(16,12))
       #  plt.subplot(2,1,1)
       #  #These y values may have to be changed
       #  plt.plot(tns, ybar_symm[0]*vel_st + vel_st, label = 'Velocity [m/s]')
       #  #plt.plot(tns,ybar_symm[1]+alpha_st, label = 'Alpha [Rad]')
       #  #plt.plot(tns,ybar_symm[2]+theta_st, label = 'Theta [Rad]')
       # # plt.plot(tns,ybar_symm[3]*vel_st/c, label = 'q [Rad/s]')
       #  plt.xlabel('Time Steps [s]')
       #  plt.ylabel(' Velocity m/s')
       #  plt.subplot(2,1,2)
       # #These y values may have to be changed
       # plt.plot(tns,deg_to_rad(flightdata['delta_e'][index_begin:index_end]), label = 'delta_e')
       # plt.xlabel('Time Steps [s]')
       # plt.ylabel('delta e [Rad]')
       # plt.show()
         
        plt.figure(figsize=(14, 10))
        plt.subplot(2,2,1)
        
        #These y values may have to be changed
        plt.plot(tns,(ybar_symm[0]*vel_st + vel_st), label = 'Simulated Response')
        plt.plot(tns, velocity[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time [s]')
        plt.ylabel('Velocity [m/s]')
        plt.subplot(2,2,2)
        #These y values may have to be changed
        plt.plot(tns,ybar_symm[1]+alpha_st, label = 'Simulated Response')
        plt.plot(tns,alpha[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time [s]')
        plt.ylabel('Alpha [Rad]')
        plt.subplot(2,2,3)
        #These y values may have to be changed
        plt.plot(tns, ybar_symm[2]+theta_st, label = 'Simulated Response')
        plt.plot(tns,theta[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time [s]')
        plt.ylabel('Theta [Rad]')
        plt.subplot(2,2,4)
        #These y values may have to be changed
        plt.plot(tns, ybar_symm[3]*vel_st/c, label = 'Simulated Response')
        plt.plot(tns, q_val[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time [s]')
        plt.ylabel('Pitch Rate [Rad/s]')
        plt.suptitle(manuever)
        plt.show()
        #This shows the error between the data provided and the simulated models
        # av_error = error(velocity[index_begin:index_end],ybar_symm[0]*vel_st + vel_st),error(alpha[index_begin:index_end],ybar_symm[1]+alpha_st), error(theta[index_begin:index_end],ybar_symm[2]+theta_st),error(q_val[index_begin:index_end],ybar_symm[3]*vel_st/c)
        # print(av_error)
        #RMS values
        RMSv = RMS(velocity[index_begin:index_end],ybar_symm[0]*vel_st + vel_st)
        RMSa = RMS(alpha[index_begin:index_end],ybar_symm[1]+alpha_st)
        RMSt = RMS(theta[index_begin:index_end],ybar_symm[2]+theta_st)
        RMSq = RMS(q_val[index_begin:index_end],ybar_symm[3]*vel_st/c)
        
        tot_error_sym.append(RMSv**2)
        tot_error_sym.append(RMSa**2)
        tot_error_sym.append(RMSt**2)
        tot_error_sym.append(RMSq**2)
        
        print('RMS velocity ', manuever, ' = ', RMSv)
        print('RMS alpha ', manuever, ' = ', RMSa)
        print('RMS theta ', manuever, ' = ', RMSt)
        print('RMS q ', manuever, ' = ', RMSq)
        print()
        print()
        
        
    else:
    #Asymm Plots
        print('-------------------ASYMMETRIC-----------------')
       # print(mataasymeig)
        plt.figure(figsize=(16,12))
        plt.subplot(2,2,1)
    # These y values may have to be changed
        plt.plot(tna,ybar_asymm[0]+yaw_beta_st, label = 'Simulated Response')
        # plt.plot(tna, yaw_beta[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time [s]')
        plt.ylabel('Beta [Rad]')
        plt.subplot(2,2,2)
         # These y values may have to be changed
        plt.plot(tna,-ybar_asymm[1] + phi_st, label = 'Simulated Response')
        plt.plot(tna,phi[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time [s]')
        plt.ylabel('Phi [Rad]')
        plt.subplot(2,2,3)
         # These y values may have to be changed
        plt.plot(tna, -(ybar_asymm[2] * (2*vel_st/b) + p_st), label = 'Simulated Response')
        plt.plot(tna,p[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time [s]')
        plt.ylabel('p Values [Rad/s]')
        plt.subplot(2,2,4)
         # These y values may have to be changed
        plt.plot(tna, -(ybar_asymm[3] * (2*vel_st/b)) + r_st, label = 'Simulated Response')
        plt.plot(tna, r[index_begin:index_end], label = 'Actual Response')
        plt.legend(loc="best")
        plt.xlabel('Time [s]')
        plt.ylabel('r Values [Rad/s]')
        plt.suptitle(manuever)
        plt.show()
        # av_error = error(phi[index_begin:index_end],ybar_asymm[1]), error(p[index_begin:index_end],-ybar_asymm[2]), error(r[index_begin:index_end],ybar_asymm[3]+r_st)
        # print(av_error)
        #RMS values
        RMSphi = RMS(phi[index_begin:index_end],ybar_asymm[1])
        RMSp = RMS(p[index_begin:index_end],-(ybar_asymm[2] * (2*vel_st/b) + p_st))
        RMSr = RMS(r[index_begin:index_end],-(ybar_asymm[3] * (2*vel_st/b) + r_st))
        
        tot_error_asym.append(RMSphi**2)
        tot_error_asym.append(RMSp**2)
        tot_error_asym.append(RMSr**2)
            
        print('RMS phi ', manuever, ' = ', RMSphi)
        print('RMS p ', manuever, ' = ', RMSp)
        print('RMS r ', manuever, ' = ', RMSr)
        print()
        print()
        
        
        #Below are the plotting functions used for when you want to plot the eigenmotions into one plot
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

tot_error_sym = np.sqrt(sum(np.asarray(tot_error_sym)))
tot_error_asym = np.sqrt(sum(np.asarray(tot_error_asym)))

print('-------------------TOTAL ERROR-------------------')
print('symmetric = ', tot_error_sym)
print('asymmetric = ', tot_error_asym)
