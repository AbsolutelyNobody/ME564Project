from math import exp

import numpy as np

# Physical Constants
KM_TO_MILES = 5/8
MILES_TO_FEET = 5280
KG_TO_LBS = 2.2
GRAVITY = 32.2 #lbf
AIR_DENSITY_GROUND = 0.002377 #slug/ft^3
AIR_DENSITY_CRUISE = 0.0012673 # slug/ft^3 at 20000 feet? Idk why they immediately assumed this
SECONDS_PER_HOUR = 60*60
FT_LBS_PER_S_PER_HP = 550

# Things that the textbook declares to be true
j = 1.15
N = 3
mu = 0.4

# Non-technical constants
CREW_MEMBER_WEIGHT = 70 # kg
CREW_PAYLOAD_WEIGHT = 10 # kg, food, water, and personal items per person

class Aircraft:
	# these are class variables, things apparently true for all aircraft
	def __init__(self, flight, airfoil):
		# Step 1: Provide a first estimate of key quantities (Ch8/pg.401)
		self.flight = flight
		self.airfoil = airfoil
		self.fuel_consumption = 2.02*10**(-7)
		self.prop_efficiency = 0.85
		self.LD_ratio = 27 # TODO: this will come from airfoild
		# look at page 406, there is reference to a calculation by Raymer
		self.airframe_fraction = 0.20 # this number is too small, was made this size to make solution exist

		cruise_ratio = 1/exp(self.fuel_consumption/self.prop_efficiency*self.flight.range/self.LD_ratio) # Eq 8.15
		# takeoff, climb, cruise, descent, landing
		# takeoff and climb can probably get closer to 1, since we have so much fuel
		self.weight_ratios = np.array([0.97, 0.985, cruise_ratio, 1, 0.995]) # fuel_fractions

		weight_f = np.prod(self.weight_ratios) # Eq 8.19

		self.fuel_fraction  = self.flight.safety * (1-weight_f) # Eq 8.8

		self.weight = self.flight.weight_stuff/(1-self.fuel_fraction-self.airframe_fraction) # Eq 8.23
		print(f"Initial weight estimate: {self.weight} lbs")

		# Sec 8.4.2
		self.wing_loading_max = self.determine_wing_loading()

		self.S = self.weight / self.wing_loading_max
		print(f"Wing Area (S): {self.S}")
		# Sec 8.4.3
		self.P, self.AR = self.determine_TWR_AR()

		print(f"Aspect Ratio: {self.AR}")


	def determine_TWR_AR(self):
		takeoff_power = self.get_power_required_takeoff()
		print(f"Power Required for Takeoff: {takeoff_power/FT_LBS_PER_S_PER_HP} hp")
		AR, climb_power = self.get_power_required_climb_and_AR()
		print(f"Power Required for Climb: {climb_power/FT_LBS_PER_S_PER_HP} hp")
		v_max_power = self.get_power_required_v_max()
		print(f"Power Required for VMax: {v_max_power/FT_LBS_PER_S_PER_HP} hp")

		max_power = np.max([takeoff_power, climb_power, v_max_power])
		print(f"Maximum of Power Req: {max_power/FT_LBS_PER_S_PER_HP} hp")
		return max_power, AR

	def get_power_required_v_max(self):
		K = 1/(4*self.airfoil.Cdo*self.LD_ratio**2)
		W_2 = self.weight * self.weight_ratios[0]* self.weight_ratios[1]
		W_mc = 0.5*W_2*(1+self.weight_ratios[2])
		T = W_mc * (0.5*AIR_DENSITY_CRUISE*self.flight.v_max**2*self.airfoil.Cdo/(W_mc/self.S)+2*K*W_mc/(AIR_DENSITY_CRUISE*self.flight.v_max**2*self.S))
		P = T*self.flight.v_max/self.prop_efficiency
		return P

	def get_power_required_climb_and_AR(self):
		K = 1/(4*self.airfoil.Cdo*self.LD_ratio**2)
		AR = 1/(np.pi*self.airfoil.e*K)
		P = self.weight/self.prop_efficiency * (self.flight.rc_max+np.sqrt(2/AIR_DENSITY_GROUND*np.sqrt(K/(3*self.airfoil.Cdo))*self.wing_loading_max)*1.155/self.LD_ratio)
		return AR, P

	def get_power_required_takeoff(self):
		# TODO: This below eq may be very incorrect, since it assumes that drag and friction are negligible, which might be wrong
		s_g_times_TWR = (1.21*(self.wing_loading_max) / (GRAVITY * AIR_DENSITY_GROUND * self.airfoil.max_cL_half_flaps)) # Eq 8.34

		v_stall_takeoff = np.sqrt(2*self.wing_loading_max/(AIR_DENSITY_GROUND*self.airfoil.max_cL_half_flaps))
		R = 6.96 * v_stall_takeoff ** 2 / GRAVITY
		take_off_angle = np.arccos(1-(self.flight.obstacle_height/R)) # Eq 6.99
		s_a = R * np.sin(take_off_angle) # Eq 8.35

		TWR_at_seventy_percent_LO = s_g_times_TWR/(self.flight.takeoff_distance-s_a) # Eq 8.36

		# Eq 8.37
		V = 0.7 * 1.1 * v_stall_takeoff
		T = TWR_at_seventy_percent_LO * self.weight
		P_required = T * V

		minimum_P = P_required / self.prop_efficiency # Eq 8.38

		return minimum_P


	def determine_wing_loading(self):
		#Wing Loading 8.4.2

		wing_loading_stall = 0.5*AIR_DENSITY_GROUND*self.flight.v_stall**2 * self.airfoil.max_cL_full_flaps # Eq 8.27
		# setting  up coefficients for quadratic formula
		a = j**2/(GRAVITY * AIR_DENSITY_GROUND * self.airfoil.max_cL_full_flaps * mu)
		b = j * N * np.sqrt(2 / (AIR_DENSITY_GROUND*self.airfoil.max_cL_full_flaps))
		c = - self.flight.landing_distance
		print(a,b,c)
		wing_loading_land  = ((-b + np.sqrt(b**2 - 4 * a * c))/(2*a))**2 # solving 8.29
		print(f"Found restrictions on wing load: {wing_loading_stall} (stall condition), {wing_loading_land} (landing_conditions)")
		driving_wing_loading = np.min([wing_loading_land, wing_loading_stall])
		print(f"Minimum W/S: {driving_wing_loading}")
		return driving_wing_loading


class Flight:
	def __init__(self, range, crew, time, safety, stall):
		self.range = range * KM_TO_MILES * MILES_TO_FEET
		self.crew = crew
		self.time = time
		# are chairs part of airframe, or are chairs in this weight
		self.weight_stuff = self.crew * (CREW_MEMBER_WEIGHT + CREW_PAYLOAD_WEIGHT) * KG_TO_LBS
		self.safety = safety
		self.v_max = 200 * KM_TO_MILES * MILES_TO_FEET / SECONDS_PER_HOUR
		self.v_stall = stall
		self.v_flare = self.v_stall * 1.23
		self.landing_distance = 2200 # ft, maybe move out the magic number
		self.takeoff_distance = 2200 # ft, maybe move out the magic number
		self.obstacle_height = 50 # ft
		self.rc_max = 1 #ft/s


class Airfoil:
	def __init__(self, filepath):
		self.airfoil_data = np.genfromtxt(filepath, delimiter=',')
		self.max_cl = self.get_max_cl()
		# using only l/L to distinguish seems bad, but they do it
		self.max_cl_full_flaps = 0.9 + self.max_cl # idk why tf this is, but -\_o_/-
		self.max_cL_full_flaps = 0.9 * self.max_cl_full_flaps # Eq 8.24

		self.max_cl_half_flaps = 0.5 + self.max_cl
		self.max_cL_half_flaps = 0.9 * self.max_cl_half_flaps

		self.e = 0.6 # TODO: I don't think this is right
		self.Cdo = 4 * 0.0043 # TODO: I don't think this is right


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
	V_STALL = 120 * KM_TO_MILES * MILES_TO_FEET / SECONDS_PER_HOUR # this is already pretty high, but we can mess with it

	# this just stores the goal, vs details of aircraft
	flight = Flight(R, CREW, FLIGHT_TIME, SAFETY_FACTOR, V_STALL)

	# running this extracts the current estimate
	airfoil = Airfoil('./xf-r1046-il-1000000.csv')

	# creating the aircraft sets up an initial guess for key quantities
	aircraft = Aircraft(flight, airfoil)

if __name__ == '__main__':
	main()