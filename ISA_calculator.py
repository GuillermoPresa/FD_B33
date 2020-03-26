import numpy as np

from math import *

# ISA Constants
R = 287.0
g0 = 9.80665
rho0 = 1.225
pressure0 = 101325.0
altn = np.array([0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0, 84852.0, 90000.0, 100000.0, 110000.0, 120000.0])
tempn = np.array([288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.95, 186.95, 201.95, 251.95])
layer_type = np.array([1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1])




#A normal layer of the atmosphere. To be used in the ISA calculator
def layer_norm(alt, i, density0, press0):
    temp0 = tempn[i]
    temp = temp0
    density = density0
    press = press0

    if alt > altn[i] and alt < altn[i + 1]:
        temp = temp0 + ((tempn[i + 1] - tempn[i]) / (altn[i + 1] - altn[i])) * (alt - altn[i])
        press = press0 * (temp / tempn[i]) ** (-g0 / (((tempn[i + 1] - tempn[i]) / (altn[i + 1] - altn[i])) * R))
        density = density0 * (temp / tempn[i]) ** ((-g0 / (((tempn[i + 1] - tempn[i]) / (altn[i + 1] - altn[i])) * R)) - 1)
    elif alt >= altn[i + 1]:
        temp = temp0 + ((tempn[i + 1] - tempn[i]) / (altn[i + 1] - altn[i])) * (altn[i + 1] - altn[i])
        press = press0 * (temp / tempn[i]) ** (-g0 / (((tempn[i + 1] - tempn[i]) / (altn[i + 1] - altn[i])) * R))
        density = density0 * (temp / tempn[i]) ** ((-g0 / (((tempn[i + 1] - tempn[i]) / (altn[i + 1] - altn[i])) * R)) - 1)
    elif alt == altn[i]:
        temp = temp0
        press = press0
        density = density0
    return density, temp, press

#An isothermal layer of the atmosphere. To be used in the ISA calculator
def iso_layer(alt, i, density0, press0):
    temp0 = tempn[i]
    temp = temp0
    density = density0
    press = press0
    if alt > altn[i] and alt < altn[i + 1]:
        press = press * exp(-g0 / (R * temp) * (alt - altn[i]))
        density = density * exp(-g0 / (R * temp) * (alt - altn[i]))
    elif alt >= altn[i + 1]:
        press = press * exp(-g0 / (R * temp) * (altn[i + 1] - altn[i]))
        density = density * exp(-g0 / (R * temp) * (altn[i + 1] - altn[i]))
    elif alt == altn[i]:
        press = press0
        density = density0
    return density, temp, press

#ISA calculator to obtain the density, temperature and pressure at a specific altitude in meters
def ISAcalc(alt):
    temp = tempn[0]
    press = pressure0
    density = rho0
    if alt < 0:
        return density, temp, press
    j = 0
    while True:

        if j >= len(layer_type):
            return density,temp,press
        init_alt = altn[j]

        try:
            top_alt = altn[j + 1]
        except IndexError:
            top_alt = altn[j] + 100
        layer_type1 = layer_type[j]

        if init_alt <= alt < top_alt:
            if layer_type1 == 1:
                density, temp, press = layer_norm(alt, j, density, press)
            elif layer_type1 == 0:
                density, temp, press = iso_layer(alt, j, density, press)

        elif alt >= top_alt:
            if layer_type1 == 1:
                density, temp, press = layer_norm(top_alt, j, density, press)
            elif layer_type1 == 0:
                density, temp, press = iso_layer(top_alt, j, press, density)
        j += 1