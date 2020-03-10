from math import *

def EAStoTAS(eas,density, density0 =1.225):
    return eas*sqrt(density0/density)

def IAStoMach(press0, density0, SeaLevelTemperature, PressureAltitude, IAS, LapseRate, gamma = 1.4):
    press = CalibratedPressure(press0, SeaLevelTemperature, PressureAltitude, LapseRate)
    return sqrt( (2/(gamma - 1)) * ((1 + (press0/press)*( (1 + (gamma - 1) / (2 * gamma) * density0/press0 * IAS**2 )**(gamma/(gamma-1)) -1 ) )**((gamma-1)/gamma) - 1 ) )

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

