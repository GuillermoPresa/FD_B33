# Citation 550 - Linear simulation
from math import *

#select original/adjusted values for validation
#if VERSION = 1 adjusted values are selected
#if VERSION = 0 original (brightspace) values are selected

VERSION = 1

# xcg = 0.25 * c

# Stationary flight condition
hp0    = 5000      	      # pressure altitude in the stationary flight condition [m]
V0     = 95            # true airspeed in the stationary flight condition [m/sec]
alpha0 = 5*pi/180            # angle of attack in the stationary flight condition [rad]
th0    = 6*pi/180            # pitch angle in the stationary flight condition [rad]

# Aircraft mass
m      = 790 + 4100*0.453592 + 9165*0.453592         # initial mass (TOW) [kg] (PAX + fuel + OEW)

# aerodynamic properties
e      =  0.8            # Oswald factor [ ]
CD0    = 0.04            # Zero lift drag coefficient [ ]
CLa    = 4.55             # Slope of CL-alpha curve [ ]

# Longitudinal stability
Cma    = -0.7            # longitudinal stabilty [ ]  FILLER
Cmde   = -1.6            # elevator effectiveness [ ] FILLER

# Aircraft geometry

S      = 30.00	          # wing area [m^2]
Sh     = 0.2 * S         # stabiliser area [m^2]
Sh_S   = Sh / S	          # [ ]
lh     = 0.71 * 5.968    # tail length [m]
c      = 2.0569	          # mean aerodynamic cord [m]
lh_c   = lh / c	          # [ ]
b      = 15.911	          # wing span [m]
bh     = 5.791	          # stabilser span [m]
A      = b ** 2 / S      # wing aspect ratio [ ]
Ah     = bh ** 2 / Sh    # stabilser aspect ratio [ ]
Vh_V   = 1	          # [ ]
ih     = -2 * pi / 180   # stabiliser angle of incidence [rad]

# Constant values concerning atmosphere and gravity

rho0   = 1.2250          # air density at sea level [kg/m^3] 
lambda1 = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15          # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)

# air density [kg/m^3]  
rho    = rho0 * pow( ((1+(lambda1 * hp0 / Temp0))), (-((g / (lambda1*R)) + 1)))   
W      = m * g            # [N]       (aircraft weight)

# Constant values concerning aircraft inertia

muc    = m / (rho * S * c)
mub    = m / (rho * S * b)
KX2    = 0.019
KZ2    = 0.042
KXZ    = 0.002
KY2    = 1.25 * 1.114

# Aerodynamic constants

Cmac   = 0                      # Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa                    # Wing normal force slope [ ]
CNha   = 2 * pi * Ah / (Ah + 2) # Stabiliser normal force slope [ ]
depsda = 4 / (A + 2)            # Downwash gradient [ ]

# Lift and drag coefficient

CL = 2 * W / (rho * V0 ** 2 * S)              # Lift coefficient [ ]
CD = CD0 + (CLa * alpha0) ** 2 / (pi * A * e) # Drag coefficient [ ]


if VERSION == 1:
    # Adjusted cit_par data for optimising the numerical model
    # Stabiblity derivatives

    CX0 = W * sin(th0) / (0.5 * rho * V0 ** 2 * S)
    CXu = -0.001  # -0.02792
    CXa = +0.2  # +0.47966		# Positive! (has been erroneously negative since 1993)
    CXadot = +0.10  # +0.08330
    CXq = -0.4  # -0.28170
    CXde = -0.03728

    CZ0 = -W * cos(th0) / (0.5 * rho * V0 ** 2 * S)
    CZu = -0.9  # -0.37616
    CZa = -5.74340
    CZadot = -0.00350
    CZq = -5.66290
    CZde = -0.69612

    Cmu = +0.13  # +0.06990
    Cmadot = +1.2  # +0.17800
    Cmq = -16  # -8.79415
    
    CYb    = -0.7500
    CYbdot =  0     
    CYp    = -0.0304
    CYr    = +0.8495
    CYda   = -0.0400
    CYdr   = +0.2300
    
    Clb    = -0.18#-0.10260
    Clp    = -0.85#-0.71085
    Clr    = +0.18#+0.23760
    Clda   = -0.23088
    Cldr   = +0.03440
    
    Cnb    =  0.1348
    Cnbdot =   0.03 #0
    Cnp    =  -0.0602
    Cnr    =  -0.2061
    Cnda   =  -0.0120
    Cndr   =  -0.0939
    
    
    cxdt = 0
    cmdt = 0
    czdt = 0
    
else:
    # original Cit_par data as downloaded from brightspace
    # Stabiblity derivatives
    
    CX0    = W * sin(th0) / (0.5 * rho * V0 ** 2 * S)
    CXu    = -0.02792
    CXa    = +0.47966		# Positive! (has been erroneously negative since 1993) 
    CXadot = +0.08330
    CXq    = -0.28170
    CXde   = -0.03728
    
    CZ0    = -W * cos(th0) / (0.5 * rho * V0 ** 2 * S)
    CZu    = -0.37616
    CZa    = -5.74340
    CZadot = -0.00350
    CZq    = -5.66290
    CZde   = -0.69612
    
    Cmu    = +0.06990
    Cmadot = +0.17800
    Cmq    = -8.79415
    
    CYb    = -0.7500
    CYbdot =  0     
    CYp    = -0.0304
    CYr    = +0.8495
    CYda   = -0.0400
    CYdr   = +0.2300
    
    Clb    = -0.10260
    Clp    = -0.71085
    Clr    = +0.23760
    Clda   = -0.23088
    Cldr   = +0.03440
    
    Cnb    =  +0.1348
    Cnbdot =   0     
    Cnp    =  -0.0602
    Cnr    =  -0.2061
    Cnda   =  -0.0120
    Cndr   =  -0.0939
    
    
    cxdt = 0
    cmdt = 0
    czdt = 0