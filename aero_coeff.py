from math import *

def IAStoTAS(eas,density, density0 =1.225):
    return eas*sqrt(density0/density)
def CAStomach(press0, press, density0, cas, gamma = 1.4):
    a = 1/(gamma-1)

    return sqrt(2*a*((1 + ((2*gamma*a)**-1)*(density0/press0)*cas**2)**((gamma*a)**-1)-1))
