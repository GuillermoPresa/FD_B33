from Cit_par import *
from math import *
import numpy as np
import matplotlib.pyplot as plt
import control
from Data_reader import *
#These matrices are in the form of (C1d/dt)*xdot + (C2d/dt)*x + C1d/dt = 0
def symmetric_matrix_c1(muc, c, vel, czadot, cmadot, ky):
    return np.array([np.array([(-2*muc*c/vel), 0, 0, 0]), np.array([0, (czadot-2*muc)*c/vel, 0, 0]), np.array([0, 0, -c/vel, 0]), np.array([0, cmadot*c/vel, 0, (-2*muc*ky**2 * c/vel)])])
def symmetric_matrix_c2(cxu,cxa,cz0,czu,cza,cx0, czq,muc, cmu, cma, cmq):
    return np.array([np.array([cxu, cxa, cz0, 0]),np.array([czu, cza, -cx0, czq+2*muc]),np.array([0,0,0,1]), np.array([cmu,cma,0,cmq])])
def symmetric_matrix_c3(cxde,cxdt,czde,czdt,cmde,cmdt):
    return np.array([np.array([cxde,cxdt]),np.array([czde,czdt]),np.array([0,0]),np.array([cmde,cmdt])])



#These matrices are for the asymmetric state space solution

def asymmetric_matrix_c1(cybdot, mub,b,vel, kx,kxz, cnbdot, kz):
    return np.array([np.array([(cybdot-2*mub)*b/vel, 0, 0, 0]),np.array([0, -b/(2*vel), 0, 0]),np.array([0,0, -4*mub*kx**2 * b/vel, 4*mub*kxz*b/vel]), np.array([cnbdot*b/vel, 0, 4*mub*kxz*b/vel, -4*mub*kz**2 * b/vel])])
def asymmetric_matrix_c2(cyb, cl, cyp, cyr, mub, clb, clp, clr, cnb, cnp, cnr):
    return np.array([np.array([cyb, cl,cyp,cyr-4*mub]),np.array([0,0,1,0]), np.array([clb, 0, clp, clr]),np.array([cnb,0,cnp,cnr])])
def asymmetric_matrix_c3(cyda, cydr, clda, cldr, cnda, cndr):
    return np.array([np.array([cyda,cydr]),np.array([0,0]),np.array([clda,cldr]),np.array([cnda,cndr])])

#this takes in the c1,c2,c3 matrices and returns the state_space of them
def abcd_solver(c1mat,c2mat,c3mat):

    a = -np.matmul(np.linalg.inv(c1mat), c2mat)
    b = -np.matmul(np.linalg.inv(c1mat), c3mat)
    c = np.identity(4)
    d = np.array([[0,0],[0,0],[0,0],[0,0]])
    return control.ss(a,b,c,d)
#this returns the ybar values from the equation ybar = C*xbar + D*ubar the output is a y, a t and an x out
def ybar(state_space,ubar,time_steps):
    return control.forced_response(state_space,time_steps,ubar)



