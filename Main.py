from tkinter import *
from matplotlib import pyplot as plt
import numpy as np
import math
import cg_pos
import Stat_mes
import aero_coeff
import ISA_calculator
from scipy import stats
import os
import csv
import time
import Data_processing
from Cit_par import A as Aspect_ratio


SeaLevelPressure = 101325 #[Pa]
SeaLevelDensity = 1.225 #[kg/m^3]
SeaLevelTemperature = 15 #[C]
SeaLevelTemperature = SeaLevelTemperature + 273.15 #[K]

#Aircraft conditions
Empty_mass = 0.453592 * 9165.0
Empty_weight = 4.44822 * 9165.0
Empty_arm = 0.0254 * 291.65
Empty_moment = 0.1129848 * 2672953.5
WingArearea = 30 #m^2
Chord = 2.0569

TotalFuelMass = 1859.729

#CREW = Passenger(mass, seat number)
Pilot1 = cg_pos.Passenger(80,1)
Pilot2 = cg_pos.Passenger(102,2)

#Group B33  = Passenger(mass, seat number)
Passenger3 = cg_pos.Passenger(69,3)
Passenger4 = cg_pos.Passenger(80,4)
Passenger5 = cg_pos.Passenger(76,5)
Passenger6 = cg_pos.Passenger(105,6)
Passenger7 = cg_pos.Passenger(95,7)
Passenger8 = cg_pos.Passenger(95,8)
Passenger10 = cg_pos.Passenger(88,10)


#Bags  = Bag(mass, baggage compartment)
Bag1 = cg_pos.Bag(0,1)
Bag2 = cg_pos.Bag(0,2)
Bag3 = cg_pos.Bag(0,2)

#Aircraft data. Also in cg_pos. Change it there if ant changes are made here
Empty_mass = 0.453592 * 9165.0
Empty_weight = 4.44822 * 9165.0
Empty_arm = 0.0254 * 291.65
Empty_moment = 0.1129848 * 2672953.5

ActualFuelMass = TotalFuelMass

                                                                                 
PassengerList = [Pilot1, Pilot2, Passenger3, Passenger4, Passenger5, Passenger6, Passenger7, Passenger8, Passenger10]
BaggageList = [Bag1, Bag2, Bag3]
PayloadList = PassengerList + BaggageList

TotalPayloadMass = 0
for Item in PayloadList:
    TotalPayloadMass = TotalPayloadMass + Item.mass




#CN = W/(0.5*density*velocity2*surface)



def Coeficients1(Static_Measurements_1, Data_reduction = False):

#     print("Executing Main: Coeficients1")
    for DataLine in Static_Measurements_1.DataLineList: #Processing1
        DataLine[7] = (ISA_calculator.ISAcalc(DataLine[0])[0])    #Append Temperature
        DataLine[8] = (ISA_calculator.ISAcalc(DataLine[0])[1])    #Append Pressure
        DataLine[9] = (ISA_calculator.ISAcalc(DataLine[0])[2])    #Append Density
        DataLine[10] = (Empty_mass +TotalPayloadMass+TotalFuelMass - DataLine[5]) #Append Mass
        print("Mass or DataLine[10] is:", DataLine[10])
        DataLine[11] = (aero_coeff.IAStoMach(SeaLevelPressure, SeaLevelDensity, SeaLevelTemperature, DataLine[0], DataLine[1]))    #Append Mach
        DataLine[12] = (DataLine[11]*aero_coeff.SpeedOfSound(DataLine[7]))    #Append TAS
        
        #if Data_reduction is True:
        #    DataLine[12] = Data_processing.red_airspeed(DataLine[0], DataLine[1], DataLine[6], (TotalFuelMass - DataLine[5]))[1]
        
        DataLine[13] = ((DataLine[10]*9.80665)/(0.5*DataLine[9]*math.pow(DataLine[12],2)*WingArearea)) #Just a placeholder until we substract the sin(a)*T of W
#         print((DataLine[10]*9.80665)/(0.5*DataLine[9]*math.pow(DataLine[12],2)*WingArearea))


    #Create input matlab for thrust.exe
    #os.remove("matlab.dat")
    PressureAltitude = np.zeros(len(Static_Measurements_1.DataLineList))
    Mach = np.zeros(len(Static_Measurements_1.DataLineList))
    FuelFlow1 = np.zeros(len(Static_Measurements_1.DataLineList))
    FuelFlow2 = np.zeros(len(Static_Measurements_1.DataLineList))
    TemperatureDifference = np.zeros(len(Static_Measurements_1.DataLineList))
    Total_array = np.zeros(5*len(Static_Measurements_1.DataLineList))

    for i in range(0,len(Static_Measurements_1.DataLineList)): #Processing1
        PressureAltitude[i] = Static_Measurements_1.DataLineList[i][0]
        Total_array[5*i] = Static_Measurements_1.DataLineList[i][0]
        Mach[i] = Static_Measurements_1.DataLineList[i][11]
        Total_array[5*i+1] = Static_Measurements_1.DataLineList[i][11]
        FuelFlow1[i] = Static_Measurements_1.DataLineList[i][3]
        Total_array[5*i+2] = Static_Measurements_1.DataLineList[i][3]
        FuelFlow2[i] = Static_Measurements_1.DataLineList[i][4]
        Total_array[5*i+3] = Static_Measurements_1.DataLineList[i][4]
        TemperatureDifference[i] = Static_Measurements_1.DataLineList[i][7] - Static_Measurements_1.DataLineList[i][6]
        Total_array[5*i+4]  = Static_Measurements_1.DataLineList[i][6] - Static_Measurements_1.DataLineList[i][7]



    dat = np.array([PressureAltitude, Mach, TemperatureDifference, FuelFlow1, FuelFlow2])
#     print("DAT array is:", dat)
    a = np.column_stack((dat))
#     print("A column stack is:", a)
    np.savetxt('matlab.dat', a, delimiter=' ')

    #Run thrust.exe

    os.startfile("C:/Users/Guille/Documents/GitHub/FD_B33/thrust.exe")

    #Read thrust.exe output
    time.sleep(0.2)
    datfile = open("thrust.dat")                                                                        
    datContent = [i.strip().split() for i in datfile.readlines()]
    for i in range(0,len(Static_Measurements_1.DataLineList)):
        Static_Measurements_1.DataLineList[i][14] = (float(datContent[i][0])+float(datContent[i][1]))
        #Static_Measurements_1.DataLineList.append(float(ThrustTupple[0])+float(ThrustTupple[1]))
    datfile.close()
    os.remove("thrust.dat")

    Cl_array = np.zeros(len(Static_Measurements_1.DataLineList))
    Cd_array = np.zeros(len(Static_Measurements_1.DataLineList))
    Alpha_array = np.zeros(len(Static_Measurements_1.DataLineList))
    i = 0


    for i in range(0,len(Static_Measurements_1.DataLineList)): #Processing2
#         print("AT i: ",i,"The AOA is: ",Static_Measurements_1.DataLineList[i][2])
#         print((Static_Measurements_1.DataLineList[i][10]*9.80665)/(0.5*Static_Measurements_1.DataLineList[i][9]*math.pow(Static_Measurements_1.DataLineList[i][12],2)*WingArearea))

        Static_Measurements_1.DataLineList[i][15] = ((math.cos(math.pi/180*Static_Measurements_1.DataLineList[i][2])*Static_Measurements_1.DataLineList[i][14])/(0.5*Static_Measurements_1.DataLineList[i][9]*math.pow(Static_Measurements_1.DataLineList[i][12],2)*WingArearea))
        Static_Measurements_1.DataLineList[i][13] = ((Static_Measurements_1.DataLineList[i][10]*9.80665-math.sin(math.pi/180*Static_Measurements_1.DataLineList[i][2])*Static_Measurements_1.DataLineList[i][14])/(0.5*Static_Measurements_1.DataLineList[i][9]*math.pow(Static_Measurements_1.DataLineList[i][12],2)*WingArearea))
        Alpha_array[i] = Static_Measurements_1.DataLineList[i][2]
        Cl_array[i] = Static_Measurements_1.DataLineList[i][13]
        Cd_array[i] = Static_Measurements_1.DataLineList[i][15]

    print(Static_Measurements_1.DataLineList)

    fig = plt.figure()
    ax1 = fig.add_subplot(411)
    ax1.set_title('Cla')
    ax1.set_xlabel('Alpha')
    ax1.set_ylabel('Cl')
    img1 = ax1.plot(Alpha_array, Cl_array)

    ax1 = fig.add_subplot(412)
    ax1.set_title('Cda')
    ax1.set_xlabel('Alpha')
    ax1.set_ylabel('Cd')
    img1 = ax1.plot(Alpha_array, Cd_array)

    ax1 = fig.add_subplot(413)
    ax1.set_title('Cl/Cd')
    ax1.set_xlabel('Cd')
    ax1.set_ylabel('Cl')
    img1 = ax1.plot(Cd_array, Cl_array)
    #ax1.set_ylim(0,1.25)
    #ax1.set_xlim(-0.000003,0.00003)    
    plt.show()


    return Alpha_array, Cl_array, Cd_array

Static_Measurements_1 = Stat_mes.DataBlock(PayloadList)

                                     #    [hp,    IAS,    a,        FFl,    FFr,    F.used,    TAT,  /*/  Temp,    Press,    Density,    Mass,    Mach,    TAS,    Cl,        Tot-Thrust,        Cd]

Static_Measurements_1.DataLineList =[   [5010,    250,    1.7,    762,    817,    286,    6.5],
                                        [6030,    220,    2.4,    660,    698,    371,    3.2],
                                        [6030,    192,    3.6,    535,    577,    406,    1.3],
                                        [6040,    161,    5.7,    454,    487,    442,    0.8]]
                                        #[6030,    132,    8.5,    417,    459,    485,    -0.5],
                                        #[6140,    117,    11.2,    436,    473,    515,    -1.5]]


for DataLine in Static_Measurements_1.DataLineList:
    for i in range(6,15):
        DataLine.append(0)
    DataLine[0] = DataLine[0]*0.3048    #ft to m
    DataLine[1] = DataLine[1]*0.514444    #Knots to m/s
    DataLine[3] = DataLine[3]*0.453592/3600    #Lb/hr to kg/s
    DataLine[4] = DataLine[4]*0.453592/3600    #Lb/hr to kg/s
    DataLine[5] = DataLine[5]*0.453592    #Lb to Kg
    DataLine[6] = DataLine[6]+273.15    #C to K


test = Coeficients1(Static_Measurements_1, True)[2]

def Coeficients2(Static_Measurements_2, Data_reduction = False):
#     print("Executing Main: Coeficients2")
    for DataLine in Static_Measurements_2.DataLineList: #Processing1
        DataLine[7+3] = (ISA_calculator.ISAcalc(DataLine[0])[0])    #Append Temperature
        DataLine[8+3] = (ISA_calculator.ISAcalc(DataLine[0])[1])    #Append Pressure
        DataLine[9+3] = (ISA_calculator.ISAcalc(DataLine[0])[2])    #Append Density
        DataLine[10+3] = (Empty_mass +TotalPayloadMass+TotalFuelMass-DataLine[8])    
        DataLine[11+3] = (aero_coeff.IAStoMach(SeaLevelPressure, SeaLevelDensity, SeaLevelTemperature, DataLine[0], DataLine[1]))    #Append Mach
        DataLine[12+3] = (DataLine[11+3]*aero_coeff.SpeedOfSound(DataLine[7+3]))    #Append TAS
        
        if Data_reduction is True:
            DataLine[15] = Data_processing.red_airspeed(DataLine[0], DataLine[1], DataLine[10], (TotalFuelMass - DataLine[8]))[1]
        
        DataLine[16] = ((DataLine[10+3]*9.80665)/(0.5*DataLine[9+3]*math.pow(DataLine[12+3],2)*WingArearea)) #Only placeholder
#         print((DataLine[10+3]*9.80665)/(0.5*DataLine[9+3]*math.pow(DataLine[12+3],2)*WingArearea)) a


    #Create input matlab for thrust.exe
    os.remove("matlab.dat")
    PressureAltitude = np.zeros(len(Static_Measurements_2.DataLineList))
    Mach = np.zeros(len(Static_Measurements_2.DataLineList))
    FuelFlow1 = np.zeros(len(Static_Measurements_2.DataLineList))
    FuelFlow2 = np.zeros(len(Static_Measurements_2.DataLineList))
    TemperatureDifference = np.zeros(len(Static_Measurements_2.DataLineList))
    Total_array = np.zeros(5*len(Static_Measurements_2.DataLineList))

    for i in range(0,len(Static_Measurements_2.DataLineList)): #Processing1
        PressureAltitude[i] = Static_Measurements_2.DataLineList[i][0]
        Total_array[5*i] = Static_Measurements_2.DataLineList[i][0]
        Mach[i] = Static_Measurements_2.DataLineList[i][11+3]
        Total_array[5*i+1] = Static_Measurements_2.DataLineList[i][11+3]
        FuelFlow1[i] = 0.048
        Total_array[5*i+2] = Static_Measurements_2.DataLineList[i][3+3]
        FuelFlow2[i] = 0.048
        Total_array[5*i+3] = Static_Measurements_2.DataLineList[i][4+3]
        TemperatureDifference[i] = Static_Measurements_2.DataLineList[i][7+3] - Static_Measurements_2.DataLineList[i][6+3]
        Total_array[5*i+4]  = Static_Measurements_2.DataLineList[i][6+3] - Static_Measurements_2.DataLineList[i][7+3]



    dat = np.array([PressureAltitude, Mach, TemperatureDifference, FuelFlow1, FuelFlow2])
#     print("DAT array is:", dat)
    a = np.column_stack((dat))
#     print("A column stack is:", a)
    np.savetxt('matlab.dat', a, delimiter=' ')

    #Run thrust.exe

    os.startfile("C:/Users/Guille/Documents/GitHub/FD_B33/thrust.exe")
    time.sleep(0.2)
    #Read thrust.exe output
    datfile = open("thrust.dat")
    datContent = [i.strip().split() for i in datfile.readlines()]
    for i in range(0,len(Static_Measurements_2.DataLineList)):
#         print(float(datContent[i][0]))
#         print("Step:", i)
        Static_Measurements_2.DataLineList[i][17] = (float(datContent[i][0])+float(datContent[i][1]))
        #Static_Measurements_2.DataLineList.append(float(ThrustTupple[0])+float(ThrustTupple[1]))
    datfile.close()
    os.remove("thrust.dat")

    Cl_array = np.zeros(len(Static_Measurements_2.DataLineList))
    Cd_array = np.zeros(len(Static_Measurements_2.DataLineList))
    Alpha_array = np.zeros(len(Static_Measurements_2.DataLineList))
    i = 0


    for i in range(0,len(Static_Measurements_2.DataLineList)): #Processing2
#         print("AT i: ",i,"The AOA is: ",Static_Measurements_2.DataLineList[i][2])
#         print((Static_Measurements_2.DataLineList[i][10+3]*9.80665)/(0.5*Static_Measurements_2.DataLineList[i][9+3]*math.pow(Static_Measurements_2.DataLineList[i][12+3],2)*WingArearea))

        Static_Measurements_2.DataLineList[i][18] = ((math.cos(math.pi/180*Static_Measurements_2.DataLineList[i][2])*Static_Measurements_2.DataLineList[i][14+3])/(0.5*Static_Measurements_2.DataLineList[i][9+3]*math.pow(Static_Measurements_2.DataLineList[i][12+3],2)*WingArearea))
        #Static_Measurements_2.DataLineList[i][13] = ((Static_Measurements_1.DataLineList[i][10+3]*9.80665-math.sin(math.pi/180*Static_Measurements_2.DataLineList[i][2])*Static_Measurements_2.DataLineList[i][14+3])/(0.5*Static_Measurements_1.DataLineList[i][9+3]*math.pow(Static_Measurements_1.DataLineList[i][12+3],2)*WingArearea))
        Alpha_array[i] = Static_Measurements_2.DataLineList[i][2]
        Cl_array[i] = Static_Measurements_2.DataLineList[i][13+3]
        Cd_array[i] = Static_Measurements_2.DataLineList[i][15+3]

    print("wight:",DataLine[10+3])

#     print(Static_Measurements_2.DataLineList)



    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax1.set_title('Cla')
    ax1.set_xlabel('Alpha')
    ax1.set_ylabel('Cl')
    img1 = ax1.plot(Alpha_array, Cl_array)

    ax1 = fig.add_subplot(312)
    ax1.set_title('Cda')
    ax1.set_xlabel('Alpha')
    ax1.set_ylabel('Cd')
    img1 = ax1.plot(Alpha_array, Cd_array)

    ax1 = fig.add_subplot(313)
    ax1.set_title('Cl/Cd')
    ax1.set_xlabel('Cd')
    ax1.set_ylabel('Cl')
    img1 = ax1.plot(Cd_array, Cl_array)
    #ax1.set_ylim(0,1.25)
    #ax1.set_xlim(-0.000003,0.00003)    
    plt.show()


    return Alpha_array, Cl_array, Cd_array

def CoeficientsCGShift(Static_Measurements_3, Xcg1, Xcg2, Data_reduction = False):
#     print("Executing Main: CoeficientsCGShift")
    for DataLine in Static_Measurements_3.DataLineList: #Processing1
        DataLine[7+3] = (ISA_calculator.ISAcalc(DataLine[0])[0])    #Append Temperature
        DataLine[8+3] = (ISA_calculator.ISAcalc(DataLine[0])[1])    #Append Pressure
        DataLine[9+3] = (ISA_calculator.ISAcalc(DataLine[0])[2])    #Append Density
        DataLine[10+3] = (Empty_mass +TotalPayloadMass+TotalFuelMass- DataLine[8])    
        DataLine[11+3] = (aero_coeff.IAStoMach(SeaLevelPressure, SeaLevelDensity, SeaLevelTemperature, DataLine[0], DataLine[1]))    #Append Mach
        DataLine[12+3] = (DataLine[11+3]*aero_coeff.SpeedOfSound(DataLine[7+3]))    #Append TAS
        
        if Data_reduction is True:
            DataLine[15] = Data_processing.red_airspeed(DataLine[0], DataLine[1], DataLine[10], (TotalFuelMass - DataLine[8]))[1]
  
        DataLine[16] = ((DataLine[10+3]*9.80665)/(0.5*DataLine[9+3]*math.pow(DataLine[12+3],2)*WingArearea))
     #   print((DataLine[10+3]*9.80665)/(0.5*DataLine[9+3]*math.pow(DataLine[12+3],2)*WingArearea))

    CN = (DataLine[13]*9.80665)/(0.5*DataLine[12]*math.pow(DataLine[15],2)*WingArearea)
    print("Xcg2-Xcg1:", Xcg2-Xcg1)
    print("Chord:", Chord)
    print("CN:", CN)
    Cmd = -1/(math.pi/180*(Static_Measurements_3.DataLineList[1][3]-Static_Measurements_3.DataLineList[0][3])) * CN * (Xcg2-Xcg1)/Chord
    
    return Cmd

                                                                                                                                                                                                                                
Static_Measurements_3 = Stat_mes.DataBlock(PayloadList)
                                    #    [hp,    IAS,    a,        de,        detr,    Fe,        FFl,    FFr,    F.used,    TAT,    Temp,    Press,    Density,    Mass,    Mach,    TAS,    Cl,        Tot-Thrust,        Cd]
Static_Measurements_3.DataLineList = [[7750,    157,    5.7,    0.3,   2.6,    1,      421,    464,    852,    -2.0],
                                      [7620,    158,    5.7,    -0.8,    2.6,    -20,    423,    467,    894,    -1.8]]

for DataLine in Static_Measurements_3.DataLineList:
    for i in range(9,18):
        DataLine.append(0)
    DataLine[0] = DataLine[0]*0.3048    #ft to m
    DataLine[1] = DataLine[1]*0.514444    #Knots to m/s
    DataLine[3+3] = DataLine[3+3]*0.453592/3600    #Lb/hr to kg/s
    DataLine[4+3] = DataLine[4+3]*0.453592/3600    #Lb/hr to kg/s
    DataLine[5+3] = DataLine[5+3]*0.453592    #Lb to Kg
    DataLine[6+3] = DataLine[6+3]+273.15    #C to K


Xcg1 = cg_pos.cg((TotalFuelMass - Static_Measurements_3.DataLineList[0][8]), PayloadList)[0]

#Moving Passenger
Passenger8_moved = cg_pos.Passenger(95,1)


PassengerList2 = [Pilot1, Pilot2, Passenger3, Passenger4, Passenger5, Passenger6, Passenger7, Passenger8_moved, Passenger10]
PayloadList2 = PassengerList2 + BaggageList

Xcg2 = cg_pos.cg((TotalFuelMass - Static_Measurements_3.DataLineList[1][8]), PayloadList2)[0]


CoeficientsCGShift(Static_Measurements_3,Xcg1,Xcg2)




                        
Static_Measurements_2 = Stat_mes.DataBlock(PayloadList)
                                    #    [hp,    IAS,    a,        de,        detr,    Fe,        FFl,    FFr,    F.used,    TAT,    Temp,    Press,    Density,    Mass,    Mach,    TAS,    Cl,        Tot-Thrust,        Cd]
Static_Measurements_2.DataLineList = [  [6910,    158,    5.7,    -0.3,   2.6,    0,    432,    470,    595,    -1.0],
                                        [7020,    149,    6.6,    -0.7,   2.6,    -7,    430,    470,    636,    -1.6],
                                        [7180,    139,    7.5,    -1.1,   2.6,    -26,    428,    467,    656,    -2.5],
                                        [7300,    129,    9.2,    -1.7,    2.6,    -42,    424,    463,    683,    -3.0],
                                        [6640,    172,    4.8,    0.2,    2.6,    29,        438,    479,    715,    0.2],
                                        [6340,    178,    4.2,    0.5,    2.6,    43,        446,    487,    736,    0.5],
                                        [5870,    188,    3.7,    0.7,    2.6,    85,        455,    496,    765,    0.8]]


for DataLine in Static_Measurements_2.DataLineList:
    for i in range(9,18):
        DataLine.append(0)
    DataLine[0] = DataLine[0]*0.3048    #ft to m
    DataLine[1] = DataLine[1]*0.514444    #Knots to m/s
    DataLine[3+3] = DataLine[3+3]*0.453592/3600    #Lb/hr to kg/s
    DataLine[4+3] = DataLine[4+3]*0.453592/3600    #Lb/hr to kg/s
    DataLine[5+3] = DataLine[5+3]*0.453592    #Lb to Kg
    DataLine[6+3] = DataLine[6+3]+273.15    #C to K


Coeficients2(Static_Measurements_2)
   
                                                                                
#-------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                DATA PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------------------------
print(Coeficients1(Static_Measurements_1))
########## STATIC MEASUREMENT 1
stat_meas1_outcomes = []

Coeficient1Output = Coeficients1(Static_Measurements_1, True)
print("Coeficient1Output", Coeficient1Output)
alpha = Coeficient1Output[0]
print("alpha:",alpha)
alpha_rad = (alpha * math.pi)/180
print("alpha_rad:",alpha_rad)
CL = Coeficient1Output[1]
print("CL:",CL)
CD = Coeficient1Output[2]
print("CD:",CD)
CL2 = np.square(CL)

# fitting lines through graphs
CL_alpha, Cl0, r_value, p_value, std_err = stats.linregress(alpha_rad, CL)
line_deriv1, Cd_0, r_value, p_value, std_err = stats.linregress(CL2[:-2], CD[:-2])
print("CLA:",CL_alpha)                                 

e = 1/(math.pi * line_deriv1 * Aspect_ratio)

#improving CL data and CD
alpha_improved = np.linspace(-0.05, 0.187, 10)
CL1 = np.linspace(-0.5, 0.5, 20)

CL_improved = CL_alpha * alpha_improved + Cl0
CD_improved = Cd_0 + (CL1**2)/(math.pi * Aspect_ratio * e)


# plotting
Data_processing.plotter_stat_meas1(alpha_improved, CL_improved, CD_improved, CL1)

# variables
stat_meas1_outcomes.append(CL_alpha)
stat_meas1_outcomes.append(Cd_0)
stat_meas1_outcomes.append(e)

print('relevant values static measurement 1: Cl_alpha, Cd_0, oswald efficiency')
print(stat_meas1_outcomes)

########## STATIC MEASURMENT 2
stat_meas2_outcomes = []

red_Ve = np.zeros(len(Static_Measurements_2.DataLineList))
print(Static_Measurements_2.DataLineList)
for i in range(0, len(Static_Measurements_2.DataLineList)):
    red_Ve[i] = Data_processing.red_airspeed(Static_Measurements_2.DataLineList[i][0], Static_Measurements_2.DataLineList[i][1], Static_Measurements_2.DataLineList[i][10], (TotalFuelMass - Static_Measurements_2.DataLineList[i][8]),Static_Measurements_2.DataLineList[i][13])[1]

el_def = np.zeros(len(Static_Measurements_2.DataLineList))
for i in range(0, len(Static_Measurements_2.DataLineList)):
    el_def[i] = math.pi/180 * (Data_processing.red_thrust(Static_Measurements_2.DataLineList[i][3], red_Ve[i], Static_Measurements_2.DataLineList[i][17], Static_Measurements_2.DataLineList[i][12]))
    
red_F_e = np.zeros(len(Static_Measurements_2.DataLineList))
for i in range(0, len(Static_Measurements_2.DataLineList)):                                                                                                                                                                                                                     
    red_F_e[i] = Data_processing.red_force(Static_Measurements_2.DataLineList[i][5], Static_Measurements_2.DataLineList[i][0], Static_Measurements_2.DataLineList[i][1], Static_Measurements_2.DataLineList[i][10], (TotalFuelMass - Static_Measurements_2.DataLineList[i][8]), Static_Measurements_2.DataLineList[i][13])

Cm_delta = CoeficientsCGShift(Static_Measurements_3, Xcg1, Xcg2, False)

alpha_deriv = Coeficients2(Static_Measurements_2, False)[0]
alpha_deriv_rad = (alpha_deriv * math.pi)/180

lin_deriv2, intersect = np.polyfit(alpha_deriv_rad, el_def, 1)
Cm_alpha = -1 * lin_deriv2 * Cm_delta


# plotting
Data_processing.plotter_stat_meas2(red_Ve, el_def, red_F_e)

# variables
stat_meas2_outcomes.append(Cm_delta)
stat_meas2_outcomes.append(Cm_alpha)


print('relevant values static measurement 1: Cl_alpha, Cd_0, oswald efficiency')
print(stat_meas1_outcomes)
print('relevant values static measurement 2: Cm_delta, Cm_alpha')
print(stat_meas2_outcomes)
                                 



                                                   
#-------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                USER INTERFACE
#-------------------------------------------------------------------------------------------------------------------------------------------------
'''
passenger_spacing_y = 40
passenger_spacing_x = 75

class MyWindow:
    def __init__(self, win):

        self.lbl_Passenger_list = Label(window, text="Passenger list", fg='red', font=("Helvetica", 16))
        self.lbl_Passenger_list.place(x=100, y=20)
    #Pilot 1
        self.LabelSeat1 = Label(window, text="Pilot 1", fg='black', font=("Helvetica", 12))
        self.LabelSeat1.place(x=50, y=75)
        self.LabelSeat1 = Label(window, text=str(Pilot1.mass) + "kg", fg='black', font=("Helvetica", 12))
        self.LabelSeat1.place(x=50+passenger_spacing_x, y=75)
    #Pilot 2
        self.LabelSeat2 = Label(window, text="Pilot 2", fg='black', font=("Helvetica", 12))
        self.LabelSeat2.place(x=50, y=75+passenger_spacing_y)
        self.WeightSeat2=Entry(window, textvariable=Pilot2.mass, bd=5)
        self.WeightSeat2.place(x=50+passenger_spacing_x, y=75+passenger_spacing_y)
    #Seat 3
        self.LabelSeat3 = Label(window, text="Seat 3", fg='black', font=("Helvetica", 12))
        self.LabelSeat3.place(x=50, y=75+2*passenger_spacing_y)
        self.WeightSeat3=Entry(window, text=Passenger3.mass, bd=5)
        self.WeightSeat3.place(x=50+passenger_spacing_x, y=75+2*passenger_spacing_y)
    #Seat 4
        self.LabelSeat4 = Label(window, text="Seat 4", fg='black', font=("Helvetica", 12))
        self.LabelSeat4.place(x=50, y=75+3*passenger_spacing_y)
        self.WeightSeat4=Entry(window, text=Passenger4.mass, bd=5)
        self.WeightSeat4.place(x=50+passenger_spacing_x, y=75+3*passenger_spacing_y)
    #Seat 5
        self.LabelSeat5 = Label(window, text="Seat 3", fg='black', font=("Helvetica", 12))
        self.LabelSeat5.place(x=50, y=75+3*passenger_spacing_y)
        self.WeightSeat5=Entry(window, text=Passenger5.mass, bd=5)
        self.WeightSeat5.place(x=50+passenger_spacing_x, y=75+3*passenger_spacing_y)
    #Seat 6
        self.LabelSeat6 = Label(window, text="Seat 4", fg='black', font=("Helvetica", 12))
        self.LabelSeat6.place(x=50, y=75+4*passenger_spacing_y)
        self.WeightSeat6=Entry(window, text=Passenger6.mass, bd=5)
        self.WeightSeat6.place(x=50+passenger_spacing_x, y=75+4*passenger_spacing_y)
    #Seat 7
        self.LabelSeat7 = Label(window, text="Seat 4", fg='black', font=("Helvetica", 12))
        self.LabelSeat7.place(x=50, y=75+5*passenger_spacing_y)
        self.WeightSeat7=Entry(window, text=Passenger7.mass, bd=5)
        self.WeightSeat7.place(x=50+passenger_spacing_x, y=75+5*passenger_spacing_y)
    #Seat 8
        self.LabelSeat8 = Label(window, text="Seat 3", fg='black', font=("Helvetica", 12))
        self.LabelSeat8.place(x=50, y=75+6*passenger_spacing_y)
        self.WeightSeat8=Entry(window, text=Passenger8.mass, bd=5)
        self.WeightSeat8.place(x=50+passenger_spacing_x, y=75+6*passenger_spacing_y)
    #Seat 9
        self.LabelSeat9 = Label(window, text="Seat 4", fg='black', font=("Helvetica", 12))
        self.LabelSeat9.place(x=50, y=75+7*passenger_spacing_y)
        self.WeightSeat9=Entry(window, text=Passenger9.mass, bd=5)
        self.WeightSeat9.place(x=50+passenger_spacing_x, y=75+7*passenger_spacing_y)
    #Seat 10
        self.LabelSeat10 = Label(window, text="Seat 4", fg='black', font=("Helvetica", 12))
        self.LabelSeat10.place(x=50, y=75+8*passenger_spacing_y)
        self.WeightSeat10=Entry(window, text=Passenger10.mass, bd=5)
        self.WeightSeat10.place(x=50+passenger_spacing_x, y=75+8*passenger_spacing_y)
    #Passenger moved
        self.lbl_Passenger_moved = Label(window, text="Passenger moved", fg='red', font=("Helvetica", 16))
        self.lbl_Passenger_moved.place(x=100, y=75+10*passenger_spacing_y)

        self.btn=Button(window, text="This is Button widget", fg='blue')
        self.btn.place(x=200, y=660)



window=Tk()
mywin=MyWindow(window)
window.title('Group B33')
window.geometry("1200x800+10+10")
window.mainloop()
'''


