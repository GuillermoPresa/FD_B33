from math import *

def IAStoTAS(eas,density, density0 =1.225):
    return eas*sqrt(density0/density)

def CAStomach(press0, press, density0, cas, gamma = 1.4):
    a = 1/(gamma-1)
    return sqrt(2*a*((1 + ((2*gamma*a)**-1)*(density0/press0)*cas**2)**((gamma*a)**-1)-1))

def SpeedOfSound(T,gamma = 1.4, R = 287.1):
    T = T + 273.15
    a = sqrt(gamma * R * T)
    return a

def StaticAirTemperature(Mach, TotalTemp, gamma = 1.4):
    return TotalTemp / (1 + ((gamma - 1)/2) * Mach)

def TAS(Mach, T):
    a = SpeedOfSound(T)
    return Mach * a

def CalibratedPressure(SeaLevelPressure, SeaLevelTemperature, PressureAltitude, LapseRate, R = 287.1, g0 = 9.80665):
	return SeaLevelPressure * (1 + LapseRate * PressureAltitude / SeaLevelTemperature)**(-g0/(LapseRate*R))

