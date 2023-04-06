from math import exp

import numpy as np

# Conversion Factors
KM_TO_MILES = 5/8
MILES_TO_FEET = 5280
KG_TO_LBS = 2.2

# Non-technical constants
CREW_MEMBER_WEIGHT = 70 # kg
CREW_PAYLOAD_WEIGHT = 10 # kg, food, water, and personal items per person

class Aircraft:
	def __init__(self, flight, airfoil):
		# Step 1: Provide a first estimate of key quantities (Ch8/pg.401)
		self.flight = flight
		self.airfoil = airfoil
		self.fuel_consumption = 2.02*10**(-7)
		self.prop_efficiency = 0.85
		self.LD_ratio = 14 # this will need to increase
		# look at page 406, there is reference to a calculation by Raymer
		self.airframe_fraction = 0.05 # this number is too small, was made this size to make solution exist
		self.planform_area = 1
		self.wetted_area = 1
		self.gravity = 32.2 #in lbf
		self.ground_air_density = 0.0765 #lb/ft^3
		self.aifoil_lift_coef = 1 # not sure how to get this yet


		cruise_ratio = 1/exp(self.fuel_consumption/self.prop_efficiency*self.flight.range/self.LD_ratio) # Eq 8.15
		# takeoff, climb, cruise, descent, landing
		# takeoff and climb can probably get closer to 1, since we have so much fuel
		self.weight_ratios = np.array([0.97, 0.985, cruise_ratio, 1, 0.995]) # fuel_fractions

		weight_f = np.prod(self.weight_ratios) # Eq 8.19

		self.fuel_fraction  = self.flight.safety * (1-weight_f) # Eq 8.8

		self.weight = self.flight.weight_stuff/(1-self.fuel_fraction-self.airframe_fraction) # Eq 8.23

		# Sec 8.4.3
		self.wing_loading = self.wetted_area / self.planform_area
		self.wing_lift_coef = 0.9*self.airfoil_lift_coef #0.9 from pg 409 as flap deflection @ 45 deg
		self.ground_roll = (1.21*(self.wing_loading) / (self.gravity * self.ground_air_density * self.wing_lift_coef *self.thrust_weight))


class Flight:
	def __init__(self, range, crew, time, safety):
		self.range = range * KM_TO_MILES * MILES_TO_FEET
		self.crew = crew
		self.time = time
		# are chairs part of airframe, or are chairs in this weight
		self.weight_stuff = self.crew * (CREW_MEMBER_WEIGHT + CREW_PAYLOAD_WEIGHT) * KG_TO_LBS
		self.safety = safety


class Airfoil:
	def __init__(self, filepath):
		self.airfoil_data = np.genfromtxt(filepath, delimiter=',')
		self.max_cl = self.get_max_cl()

	def get_max_cl(self):
		max_cl = 0
		for line in self.airfoil_data:
			cl = line[1]
			if cl>max_cl:
				max_cl = cl
		return max_cl

def main():
	# Flight information
	R = 40000 # km
	CREW = 2 # people
	FLIGHT_TIME = 10 # days
	SAFETY_FACTOR = 1.01

	# this just stores the goal, vs details of aircraft
	flight = Flight(R, CREW, FLIGHT_TIME, SAFETY_FACTOR)

	# running this extracts the current estimate
	airfoil = Airfoil('./xf-r1046-il-1000000.csv')

	# creating the aircraft sets up an initial guess for key quantities
	aircraft = Aircraft(flight, airfoil)
	print(f"Initial Weight Estimate: {aircraft.weight} lbs")


if __name__ == '__main__':
	main()