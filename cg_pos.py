import math

'''
Fix datum lines.
'''

#This module computes the cg and returns it as a funct of fuel consumed

seat_positions_x = [3.3274, 3.3274, 5.4356, 5.4356, 6.3754 , 6.3754, 7.3152, 7.3152, 3.404 ,4.318] #Its in meters
seat_positions_y = [0,0,0,0,0,0,0,0,0,0] #Its in meters

baggage_compartment_x =	[1.8796, 8.1534, 8.5852]
baggage_compartment_y = [0,0,0]

class Passenger:
	
	def __init__(self, mass, seat):
		self.mass = mass																					
		self.seat = seat
		self.position_x = seat_positions_x[seat-1]
		self.position_y = seat_positions_y[seat-1]

class Bag:
	
	def __init__(self, mass, baggage_compartment):
		self.mass = mass
		self.seat = baggage_compartment
		self.position_x = baggage_compartment_x[baggage_compartment-1]
		self.position_y = baggage_compartment_y[baggage_compartment-1]

def payload_cg(payload_list):
	total_mass = 0
	total_moment_x = 0
	total_moment_y = 0

	for item in payload_list:
		total_mass = total_mass + item.mass
		total_moment_x = total_moment_x + item.position_x * item.mass
		total_moment_y = total_moment_y + item.position_y * item.mass

	x_pos = total_moment_x/total_mass 
	y_pos = total_moment_y/total_mass 

	return total_mass, x_pos, y_pos
	
#In inch-pound/100
Fuel_moment_list = [0, 298.16, 591.18, 879.08, 1165.42, 1448.40, 1732.53, 2014.80, 2298.84, 2581.92, 2866.30, 3150.18, 3434.52, 3718.52, 4003.23, 4287.76, 4572.24, 4856.56, 5141.16, 5425.64, 5709.90, 5994.04, 6278.47, 6562.82, 6846.96, 7131.00, 7415.33, 7699.60, 7984.34, 8269.06, 8554.05, 8839.04, 9124.80, 9410.62, 9696.97, 9983.40, 10270.08, 10556.84, 10843.87, 11131.00, 11418.20, 11705.50, 11993.31, 12281.18, 12569.04, 12856.86, 13144.73, 13432.48, 13720.56, 14008.46, 14320.34]
#In pounds
Fuel_mass_list = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5008]

#Input actual fuel capacity in kg
def Fuel_cg(Actual_fuel_mass):
	if Actual_fuel_mass < 0:
		raise ValueError("executing Fuel_cg(). Actual fuel mass is negative")
	#In inch and pounds
	Fuel_index = 2.20462*Actual_fuel_mass/100
	Actual_fuel_moment = 100*(Fuel_moment_list[int(math.trunc(Fuel_index))] + (Fuel_index - math.trunc(Fuel_index)) * (Fuel_moment_list[int(math.trunc(Fuel_index))+1] - Fuel_moment_list[int(math.trunc(Fuel_index))]))
	#In nm:
	Actual_fuel_moment = 0.1129848*Actual_fuel_moment
	Fuel_cg_x = Actual_fuel_moment/Actual_fuel_mass/9.80665
	return 	Fuel_cg_x, Actual_fuel_moment

def cg(Actual_fuel_mass, payload_list): #In kg
	Passengers_mass, Passengers_cg_x, Passengers_cg_y = payload_cg(payload_list)
	Zero_fuel_mass = Empty_mass + Passengers_mass
	Zero_fuel_cg_x = (Empty_mass*Empty_arm + Passengers_cg_x*Passengers_mass) / Zero_fuel_mass
	Fuel_cg_x = 0
	Fuel_cg_x, Actual_fuel_moment = Fuel_cg(Actual_fuel_mass)
	Total_mass = Zero_fuel_mass + Actual_fuel_mass
	Total_cg_x = (Fuel_cg_x*Actual_fuel_mass + Zero_fuel_cg_x*Zero_fuel_mass) / Total_mass
	#print("The total mass is:", Total_mass, "The total x_cg is:", Total_cg_x, "or (in inch):", Total_cg_x*39.3701)
	return 	Total_cg_x, Total_mass

#CREW = Passenger(mass, seat number)
Pilot1 = Passenger(61,1)
Pilot2 = Passenger(61,2)
Coordinator = Passenger(61,10)

#Group B33  = Passenger(mass, seat number)
Guille = Passenger(61,3)
Igor = Passenger(61,4)
Bob = Passenger(61,5)
Daniel = Passenger(61,6)
Luke = Passenger(61,7)
Flip = Passenger(61,8)

#Bags  = Bag(mass, baggage compartment)
Bag1 = Bag(5,1)
Bag2 = Bag(5,2)
Bag3 = Bag(5,2)

#Aircraft data
Empty_mass = 0.453592 * 9165.0
Empty_weight = 4.44822 * 9165.0
Empty_arm = 0.0254 * 291.65
Empty_moment = 0.1129848 * 2672953.5

if abs(Empty_moment-Empty_weight*Empty_arm) > 50:
	raise ValueError("Aircraft cg readings make no sense. abs(Moment-Mass*Arm) > 50")


passenger_list = [Pilot1, Pilot2, Coordinator, Guille, Igor, Bob, Daniel, Luke, Flip]
Baggage_list = [Bag1, Bag2, Bag3]
payload_list = passenger_list + Baggage_list
Passengers_mass, Passengers_cg_x, Passengers_cg_y = payload_cg(payload_list)

Zero_fuel_mass = Empty_mass + Passengers_mass
Zero_fuel_cg_x = (Empty_mass*Empty_arm + Passengers_cg_x*Passengers_mass) / Zero_fuel_mass
#print("The Zero fuel mass is:", Zero_fuel_mass, "The zero fuel x_cg is:", Zero_fuel_cg_x, "or (in inch):", Zero_fuel_cg_x*39.3701 )

cg(100/2.2, payload_list)

