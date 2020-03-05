#This module computes the cg and returns it as a funct of fuel consumed

seat_positions_x = [0,0,0,0,0,0,0,0,0,0]
seat_positions_y = [0,0,0,0,0,0,0,0,0,0]


class Passenger:
	
	def __init__(self, mass, seat):
		self.mass = mass
		self.seat = seat
		self.seat_position_x = seat_positions_x[seat-1]
		self.seat_position_y = seat_positions_y[seat-1]

def passengers_cg(passenger_list):
	total_mass = 0
	total_moment_x = 0
	total_moment_y = 0

	for passenger in passenger_list:
		total_mass = total_mass + passenger.mass
		total_moment_x = total_moment_x + passenger.seat_position_x * passenger.mass
		total_moment_y = total_moment_y + passenger.seat_position_y * passenger.mass

	x_pos = total_moment_x/total_mass 
	y_pos = total_moment_y/total_mass 

	return total_mass, x_pos, y_pos
	

#CREW
Pilot1 = Passenger(61,1)
Pilot2 = Passenger(61,2)
Coordinator = Passenger(61,10)

#Group B33
Guille = Passenger(61,3)
Igor = Passenger(61,3)
Bob = Passenger(61,3)
Daniel = Passenger(61,3)
Luke = Passenger(61,3)
Flip = Passenger(61,3)

passenger_list = [Pilot1, Pilot2, Coordinator, Guille, Igor, Bob, Daniel, Luke, Flip]
Passengers_mass, Passengers_cg_x, Passengers_cg_y = passengers_cg(passenger_list)

#Aircraft cg readings

Empty_mass =  9165.0
Empty_arm = 291.65
Empty_moment = 2672953.5

if abs(Empty_moment-Empty_mass*Empty_arm) > 50:
	raise ValueError("Aircraft cg readings make no sense. abs(Moment-Mass*Arm) > 50")

Zero_fuel_mass = Empty_mass + Passengers_mass
Zero_fuel_cg_x = (Empty_moment + Passengers_cg_x*Passengers_mass) / Zero_fuel_mass
print("The Zero fuel mass is:", Zero_fuel_mass, "The zero fuel x_cg is:", Zero_fuel_cg_x )



