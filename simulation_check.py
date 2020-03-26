from Cit_par import *
from Data_reader import *
import matplotlib.pyplot as plt
from math import *
from ISA_calculator import ISAcalc
from state_space import *
vel_st = V0
theta_st = th0
density_st = rho0
W_new = m
alpha_st = alpha0
CL_new = 2 * W_new / (density_st * vel_st ** 2 * S)
CD_new = CD0 + (CLa * alpha_st ** 2) / (pi * A * e)
CX0_new = W_new * sin(theta_st) / (0.5 * density_st * vel_st ** 2 * S)
CZ0_new = -W_new * cos(theta_st) / (0.5 * density_st * vel_st ** 2 * S)
symmatc1 = symmetric_matrix_c1(muc, c, vel_st, CZadot, Cmadot, KY2)
symmatc2 = symmetric_matrix_c2(CXu, CXa, CZ0_new, CZu, CZa, CX0, CZq, muc, Cmu, Cma, Cmq)
symmatc3 = symmetric_matrix_c3(CXde, cxdt, CZde, czdt, Cmde, cmdt)
asymmatc1 = asymmetric_matrix_c1(CYbdot, mub, b, vel_st, KX2, KXZ, Cnbdot, KZ2)
asymmatc2 = asymmetric_matrix_c2(CYb, CL_new, CYp, CYr, mub, Clb, Clp, Clr, Cnb, Cnp, Cnr)
asymmatc3 = asymmetric_matrix_c3(CYda, CYdr, Clda, Cldr, Cnda, Cndr)
ybar_symm,time1 = ybar(symabcd_solver(symmatc1,symmatc2,symmatc3),0,np.linspace(0,50,300),[0,alpha0,0,0])
ybar_asymm,time2 = ybar(asymabcd_solver(asymmatc1,asymmatc2,asymmatc3),0,np.linspace(0,50,300),[0,alpha0,alpha0,0])
# yaw_beta = [alpha0]
# for i in range(len(r)-1):
#     yaw_beta.append(yaw_beta[i] + r[i])

# plt.figure(figsize=(14, 10))
# plt.subplot(2, 2, 1)
# plt.plot(time2, (ybar_asymm[:, 0]))
# plt.xlabel('Time [s]')
# plt.ylabel('Dimensionless Beta')
# plt.subplot(2, 2, 2)
# plt.plot(time2, ybar_asymm[:,1])
# plt.xlabel('Time [s]')
# plt.ylabel('Dimensionless Phi')
# plt.subplot(2,2,3)
# plt.plot(time2, ybar_asymm[:,2])
# plt.xlabel('Time [s]')
# plt.ylabel('Dimensionless roll Rate')
# plt.subplot(2,2,4)
# plt.plot(time2, ybar_asymm[:,3])
# plt.xlabel('Time [s]')
# plt.ylabel('Dimensionless Yaw Rate')
# plt.subplot(2, 2, 2)
# plt.plot(tns, ybar_symm[0:, 1] + alpha_st, label='Simulated Response')
# plt.plot(tns, alpha[index_begin:index_end], label='Actual Response')
# plt.legend(loc="best")
# plt.xlabel('Time Steps')
# plt.ylabel('Alpha [Rad]')
# plt.subplot(2, 2, 3)
# plt.plot(tns, ybar_symm[0:, 2] + theta_st, label='Simulated Response')
# plt.plot(tns, theta[index_begin:index_end], label='Actual Response')
# plt.legend(loc="best")
# plt.xlabel('Time Steps')
# plt.ylabel('Theta [Rad]')
# plt.subplot(2, 2, 4)
# plt.plot(tns, ybar_symm[0:, 3] * (vel_st / c), label='Simulated Response')
# plt.plot(tns, q_val[index_begin:index_end], label='Actual Response')
# plt.legend(loc="best")
# plt.xlabel('Time Steps')
# plt.ylabel('Pitch Rate [Rad/s]')
# plt.suptitle(manuever)
plt.show()
