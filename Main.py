from tkinter import *
from matplotlib import pyplot as plt
import numpy as np
import math
import cg_pos
import Stat_mes
import aero_coeff
import ISA_calculator
import scipy.io
import os
import csv





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

TotalFuelMass = 1000

#CREW = Passenger(mass, seat number)
Pilot1 = cg_pos.Passenger(61,1)
Pilot2 = cg_pos.Passenger(61,2)

#Group B33  = Passenger(mass, seat number)
Passenger3 = cg_pos.Passenger(61,3)
Passenger4 = cg_pos.Passenger(61,4)
Passenger5 = cg_pos.Passenger(61,5)
Passenger6 = cg_pos.Passenger(61,6)
Passenger7 = cg_pos.Passenger(61,7)
Passenger8 = cg_pos.Passenger(61,8)
Passenger9 = cg_pos.Passenger(61,9)
Passenger10 = cg_pos.Passenger(61,10)


#Bags  = Bag(mass, baggage compartment)
Bag1 = cg_pos.Bag(5,1)
Bag2 = cg_pos.Bag(5,2)
Bag3 = cg_pos.Bag(5,2)

#Aircraft data
Empty_mass = 0.453592 * 9165.0
Empty_weight = 4.44822 * 9165.0
Empty_arm = 0.0254 * 291.65
Empty_moment = 0.1129848 * 2672953.5

ActualFuelMass = 1000

																				 
PassengerList = [Pilot1, Pilot2, Passenger3, Passenger4, Passenger5, Passenger6, Passenger7, Passenger8, Passenger9, Passenger10]
BaggageList = [Bag1, Bag2, Bag3]
PayloadList = PassengerList + BaggageList

TotalPayloadMass = 0
for Item in PayloadList:
	TotalPayloadMass = TotalPayloadMass + Item.mass

Xcg1 = cg_pos.cg(ActualFuelMass, PayloadList)

#Moving Passenger
Passenger3 = cg_pos.Passenger(1,3)
Passenger4 = cg_pos.Passenger(61,5)


PassengerList = [Pilot1, Pilot2, Passenger3, Passenger4, Passenger5, Passenger6, Passenger7, Passenger8, Passenger9, Passenger10]
BaggageList = [Bag1, Bag2, Bag3]
PayloadList = PassengerList + BaggageList

Xcg2 = cg_pos.cg(ActualFuelMass, PayloadList)


#CN = W/(0.5*density*velocity2*surface)

Static_Measurements_1 = Stat_mes.DataBlock(PayloadList)

									#	[hp,	IAS,	a,		FFl,	FFr,	F.used,	TAT,	Temp,	Press,	Density,	Mass,	Mach,	TAS,	Cl,		Tot-Thrust,		Cd]

Static_Measurements_1.DataLineList =[	[5010,	249,	1.7,	798,	813,	360,	12.5],
										[5020,	221,	2.4,	673,	682,	412,	10.5],
										[5020,	192,	3.6,	561,	579,	447,	8.8],
										[5030,	163,	5.4,	463,	484,	478,	7.2],
										[5020,	130,	8.7,	443,	467,	532,	6],
										[5110,	118,	10.6,	474,	499,	570,	5.2]]


for DataLine in Static_Measurements_1.DataLineList: #Processing1
	DataLine[0] = DataLine[0]*0.3048	#ft to m
	DataLine[1] = DataLine[1]*0.514444	#Knots to m/s
	DataLine[3] = DataLine[3]*0.453592/3600	#Lb/hr to kg/s
	DataLine[4] = DataLine[4]*0.453592/3600	#Lb/hr to kg/s
	DataLine[5] = DataLine[5]*0.453592	#Lb to Kg
	DataLine[6] = DataLine[6]+273.15	#C to K
	DataLine.append(ISA_calculator.ISAcalc(DataLine[0])[0])	#Append Temperature
	DataLine.append(ISA_calculator.ISAcalc(DataLine[0])[1])	#Append Pressure
	DataLine.append(ISA_calculator.ISAcalc(DataLine[0])[2])	#Append Density
	#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DataLine.append(Empty_mass +TotalPayloadMass+TotalFuelMass)	#Append Mass. MISSING SUBSTRACT CONSUMED FUEL MASS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DataLine.append(aero_coeff.IAStoMach(SeaLevelPressure, SeaLevelDensity, SeaLevelTemperature, DataLine[0], DataLine[1]))	#Append Mach
	DataLine.append(DataLine[11]*aero_coeff.SpeedOfSound(DataLine[7]))	#Append TAS
	DataLine.append((DataLine[10]*9.80665)/(0.5*DataLine[9]*math.pow(DataLine[12],2)*WingArearea))
	print((DataLine[10]*9.80665)/(0.5*DataLine[9]*math.pow(DataLine[12],2)*WingArearea))


#Create input matlab for thrust.exe
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
print("DAT array is:", dat)
a = np.column_stack((dat))
#print("A column stack is:", a)
np.savetxt('matlab.dat', a, delimiter=' ')

#Run thrust.exe

os.startfile("C:/Users/Guille/Documents/GitHub/FD_B33/thrust.exe")

#Read thrust.exe output

datContent = [i.strip().split() for i in open("thrust.dat").readlines()]
for i in range(0,len(Static_Measurements_1.DataLineList)):
	Static_Measurements_1.DataLineList[i].append(float(datContent[i][0])+float(datContent[i][1]))
	#Static_Measurements_1.DataLineList.append(float(ThrustTupple[0])+float(ThrustTupple[1]))

Cl_array = np.zeros(len(Static_Measurements_1.DataLineList))
Cd_array = np.zeros(len(Static_Measurements_1.DataLineList))
Alpha_array = np.zeros(len(Static_Measurements_1.DataLineList))
i = 0


for i in range(0,len(Static_Measurements_1.DataLineList)): #Processing2
	print("AT i: ",i,"The AOA is: ",Static_Measurements_1.DataLineList[i][2])
	print((Static_Measurements_1.DataLineList[i][10]*9.80665)/(0.5*Static_Measurements_1.DataLineList[i][9]*math.pow(Static_Measurements_1.DataLineList[i][12],2)*WingArearea))

	Static_Measurements_1.DataLineList[i].append((Static_Measurements_1.DataLineList[i][13])/(0.5*Static_Measurements_1.DataLineList[i][9]*math.pow(Static_Measurements_1.DataLineList[i][12],2)*WingArearea))
	Alpha_array[i] = Static_Measurements_1.DataLineList[i][2]
	Cl_array[i] = Static_Measurements_1.DataLineList[i][13]
	Cd_array[i] = Static_Measurements_1.DataLineList[i][15]



print(Static_Measurements_1.DataLineList)



fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title('Cla')
ax1.set_xlabel('Alpha')
ax1.set_ylabel('Cl')
img1 = ax1.plot(Cl_array, Alpha_array)

ax1 = fig.add_subplot(122)
ax1.set_title('Cda')
ax1.set_xlabel('Alpha')
ax1.set_ylabel('Cd')
img1 = ax1.plot(Cd_array, Alpha_array)

ax1 = fig.add_subplot(133)
ax1.set_title('Cl/Cd')
ax1.set_xlabel('Cl')
ax1.set_ylabel('Cd')
img1 = ax1.plot(Cd_array, Cl_array)
ax1.set_ylim(0,1.25)
ax1.set_xlim(-0.000003,0.00003)	
plt.show()
	

						
Static_Measurements_ETC = Stat_mes.DataBlock(PayloadList)

Static_Measurements_ETC.DataLineList = [[6060,	161,	5.3,	0,		2.8,	0,		462,	486,	664,	5.5],
										[6350,	150,	6.3,	-0.4,	2.8,	-23,	458,	482,	694,	4.5],
										[6550,	140,	7.3,	-0.9,	2.8,	-29,	454,	477,	730,	3.5],
										[6880,	130,	8.5,	-1.5,	2.8,	-46,	449,	473,	755,	2.5],
										[6160,	173,	4.5,	0.4,	2.8,	26,		465,	489,	798,	5.0],
										[5810,	179,	4.1,	0.6,	2.8,	40,		472,	496,	825,	6.2],
										[5310,	192,	3.4,	1,		2.8,	83,		482,	505,	846,	8.2]]

	
																						 
AvgDataLine_ETC = [0] * 10
for DataLine in Static_Measurements_ETC.DataLineList:
	for i in range(0,len(DataLine)):
		AvgDataLine_ETC[i] = AvgDataLine_ETC[i]+DataLine[i]/len(Static_Measurements_ETC.DataLineList)







																					
#-------------------------------------------------------------------------------------------------------------------------------------------------
#																USER INTERFACE
#-------------------------------------------------------------------------------------------------------------------------------------------------

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



