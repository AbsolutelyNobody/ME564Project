from math import exp
import json

from matplotlib import pyplot as plt
import numpy as np

# on page 448  of page 458-398 (50/60 pages done)
# Physical Constants
KM_TO_MILES = 5/8
MILES_TO_FEET = 5280
KG_TO_LBS = 2.2
GRAVITY = 32.2 #lbf
AIR_DENSITY_GROUND = 0.002377 #slug/ft^3
AIR_DENSITY_CRUISE = 0.0013553 # slug/ft^3 at 18000 feet
SECONDS_PER_HOUR = 60*60
FT_LBS_PER_S_PER_HP = 550
AV_GAS_DENSITY = 5.64 # lb / gal
CUBIC_FEET_PER_GALLON = 0.134 # ft^3/gal
INCHES_PER_FOOT = 12

# Things that the textbook declares to be true
j = 1.15
N = 3
mu_brakes = 0.4
mu_no_brakes = 0.04

# Non-technical constants
CREW_MEMBER_WEIGHT = 70 # kg
CREW_PAYLOAD_WEIGHT = 10 # kg, food, water, and personal items per person

class Aircraft:
	# these are class variables, things apparently true for all aircraft
	def __init__(self, flight, airfoil, engine):
		# Step 1: Provide a first estimate of key quantities (Ch8/pg.401)
		self.flight = flight
		self.airfoil = airfoil
		self.engine = engine
		self.fuel_consumption = 2.02*10**(-7)
		self.prop_efficiency = 0.85
		self.fuselage_diam = 4.5
		self.LD_ratio = 26 # TODO: this will come from airfoild
		self.taper_ratio = 0.3
		# look at page 406, there is reference to a calculation by Raymer
		self.airframe_fraction = 0.25 # this number is too small, was made this size to make solution exist

		cruise_ratio = 1/exp(self.fuel_consumption/self.prop_efficiency*self.flight.range/self.LD_ratio) # Eq 8.15
		# takeoff, climb, cruise, descent, landing
		# takeoff and climb can probably get closer to 1, since we have so much fuel
		self.weight_ratios = np.array([0.97, 0.985, cruise_ratio, 1, 0.995]) # fuel_fractions

		weight_f = np.prod(self.weight_ratios) # Eq 8.19

		self.fuel_fraction  = self.flight.safety * (1-weight_f) # Eq 8.8

		self.weight = self.flight.weight_stuff/(1-self.fuel_fraction-self.airframe_fraction) # Eq 8.23
		print(f"Initial weight estimate: {self.weight} lbs")
		self.fuel_weight_est = self.weight * self.fuel_fraction
		print(f"Initial fuel mass estimate: {self.fuel_weight_est} lbs")

		print(f"Initial fuel vol estimate: {self.weight * self.fuel_fraction/AV_GAS_DENSITY} gal")
		self.fuel_vol_est = self.weight * self.fuel_fraction/AV_GAS_DENSITY * CUBIC_FEET_PER_GALLON
		print(f"Initial fuel vol estimate: {self.fuel_vol_est} ft^3")

		# Sec 8.4.2
		self.wing_loading_max = self.determine_wing_loading()

		self.S = self.weight / self.wing_loading_max
		print(f"Wing Area (S): {self.S}")


		# Sec 8.4.3
		self.P, self.AR = self.determine_TWR_AR()

		print(f"Aspect Ratio: {self.AR}")

		self.b, self.cr, self.ct, self.y_bar, self.c_bar = self.determine_wing_geometry(self.S, self.AR, self.taper_ratio)

		print("\nSingle wing geometry")
		print(f"semi-span: {self.b/2}")
		print(f"root chord: {self.cr}")
		print(f"tip chord: {self.ct}")
		print(f"y bar: {self.y_bar}")
		print(f"c bar: {self.c_bar}")

		self.wing_location = self.select_wing_location()
		print(f"wing location: {self.wing_location}")

		# checkpoint: 8.6.4 - Fuselage Configuration
		wing_fuel = self.fuel_in_wings()
		print(f"\nfuel storeable in wings: {wing_fuel} ft^3 or {wing_fuel/self.fuel_vol_est}")
		print("no fuel in wings, this isn't significant")

		cog_no_wing, self.cog_wing, self.body_lengths, wing_center_geo = self.place_components()
		print(f"cog without wing included: {cog_no_wing} ft")
		print(f"cog with wing included: {self.cog_wing} ft")
		print(f"front of engine to back of fuel tank: {self.body_lengths[4]} ft")
		print(f"front of engine to middle of wing: {wing_center_geo} ft")

		self.Sht, self.Svt = self.calculate_tail()

		# next step is propellor, 8.6.6
		self.prop_diam = self.propellor_calc()

		self.landing_gear_and_wing_placement_calc()

		self.other_configuration_stuff()

		self.weight = self.updated_weight_estimate()

		self.updated_performance()

# Calculation Functions Below

# Sec 8.4.1
class Airfoil:
	def __init__(self, performance_file, geometry_file):
		self.airfoil_performance = np.genfromtxt(performance_file, delimiter=',')
		self.airfoil_geometry = np.genfromtxt(geometry_file, delimiter=',')
		self.max_cl = self.get_max_cl()
		self.max_cL = self.max_cl * 0.9
		# using only l/L to distinguish seems bad, but they do it
		self.max_cl_full_flaps = 0.9 + self.max_cl # idk why tf this is, but -\_o_/-
		self.max_cL_full_flaps = 0.9 * self.max_cl_full_flaps # Eq 8.24

		self.max_cl_half_flaps = 0.5 + self.max_cl
		self.max_cL_half_flaps = 0.9 * self.max_cl_half_flaps

		self.e = 0.6 # TODO: I don't think this is right
		self.Cdo = 4 * 0.0043 # TODO: I don't think this is right


	def get_max_cl(self):
		max_cl = 0
		for line in self.airfoil_performance:
			cl = line[1]
			if cl>max_cl:
				max_cl = cl
		return max_cl

	def plot_airfoil(self):
		plt.scatter(self.airfoil_geometry[:,0], self.airfoil_geometry[:,1])
		plt.scatter(self.airfoil_geometry[:,0], self.airfoil_geometry[:,2])
		plt.show()

# Sec 8.4.2
	def determine_wing_loading(self):

		wing_loading_stall = 0.5 * AIR_DENSITY_GROUND * self.flight.v_stall**2 * self.airfoil.max_cL_full_flaps # Eq 8.27
		# setting  up coefficients for the quadratic formula
		a = j**2/(GRAVITY * AIR_DENSITY_GROUND * self.airfoil.max_cL_full_flaps * mu_brakes)
		b = j * N * np.sqrt(2 / (AIR_DENSITY_GROUND * self.airfoil.max_cL_full_flaps))
		c = - self.flight.landing_distance

		wing_loading_land  = ((-b + np.sqrt(b**2 - 4 * a * c))/(2*a))**2 # solving 8.29
		print(f"Found restrictions on wing load: {wing_loading_stall} (stall condition), {wing_loading_land} (landing_conditions)")
		
		driving_wing_loading = np.min([wing_loading_land, wing_loading_stall])
		print(f"Minimum W/S: {driving_wing_loading}")
		
		return driving_wing_loading

		# Sec 8.4.3
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
		print(K)
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
		v_stall_takeoff = np.sqrt(2*self.wing_loading_max/(AIR_DENSITY_GROUND*self.airfoil.max_cL_half_flaps))
		R = 6.96 * v_stall_takeoff ** 2 / GRAVITY
		take_off_angle = np.arccos(1-(self.flight.obstacle_height/R)) # Eq 6.99
		s_a_profile = R * np.sin(take_off_angle) # Eq 8.35
		s_a_climb = v_stall_takeoff*(self.flight.obstacle_height/self.flight.rc_max)
		print(f"s_a_profile: {s_a_profile}")
		print(f"s_a_climb: {s_a_climb}")
		s_a = np.max([s_a_profile, s_a_climb])
		print(f"R:{R}")
		print(f"angle:{take_off_angle*180/np.pi}")

		# TWR_at_seventy_percent_LO = s_g_times_TWR/(self.flight.takeoff_distance-s_a) # Eq 8.36
		s_g = self.flight.takeoff_distance - s_a
		a = 0.5 * AIR_DENSITY_GROUND * v_stall_takeoff ** 2 * self.S
		L = self.airfoil.max_cL_half_flaps * a
		D = L / self.LD_ratio
		TWR_at_seventy_percent_LO = (D/self.weight)+(mu_no_brakes*(1-L/self.weight))+(1.21*self.weight/self.S)/((s_g-1.1*N*np.sqrt((2*self.weight)/(AIR_DENSITY_GROUND*self.S*self.airfoil.max_cL_half_flaps)))*GRAVITY*AIR_DENSITY_GROUND*self.airfoil.max_cL_half_flaps)
		# Eq 8.37
		V = 0.7 * 1.1 * v_stall_takeoff
		T = TWR_at_seventy_percent_LO * self.weight
		P_required = T * V

		minimum_P = P_required / self.prop_efficiency # Eq 8.38

		return minimum_P

# Sec 8.5 
	def updated_performance(self):
		wing_loading = self.weight / self.S
		power_loading = self.weight / self.engine.nominal_power

		print("\n\nFinal Evaluation of Perfomance")

		print("\nKey Performance numbers")
		print(f"Wing Loading {wing_loading} lb/ft^2")
		print(f"Power Loading {power_loading} lb/hp")

		K = 1/(4*self.airfoil.Cdo*self.LD_ratio**2)
		print("\nKey Aero Parameters")
		print(f"Cdo {self.airfoil.Cdo}")
		print(f"K {K}")
		print(f"Cl_max {self.airfoil.max_cL}")
		print(f"Cl_max_flaps {self.airfoil.max_cL_full_flaps}")

		print("\nOther Paramaters")
		# solve nightmare shit
		e = 2*K*self.S/AIR_DENSITY_CRUISE*(self.weight/self.S)**2
		d = -self.engine.nominal_power
		c = 0
		b = 0
		a = 0.5*AIR_DENSITY_CRUISE*self.airfoil.Cdo*self.S
		v_max = np.polynomial.polynomial.polyroots([e,d,c,b,a])
		print(f"V_max: {v_max}")
	
	    # Sec 8.6.2 justification in the report
	def select_wing_location(self):
		# I already wrote the call and print bit in init before I realized there is no math
		return "mid-wing"

	def determine_wing_geometry(self, S, AR, taper_ratio):
		b = np.sqrt(S * AR)
		cr = 2 * S / ((taper_ratio + 1) * b) #root chord length
		ct = taper_ratio * cr # tip chord length

		y_bar = (b / 6) * ((1 + 2 * taper_ratio) / (1 + taper_ratio))
		c_bar = (2 / 3) * cr * ((1 + taper_ratio + taper_ratio ** 2) / (1 + taper_ratio))
		return b, cr, ct, y_bar, c_bar

# Sec 8.6.3 & 8.6.4
	def place_components(self):
		# assuming rectangular for now, front (engine) section, with datum point at the nose
		engine_front = 0
		engine_back = engine_front + self.engine.dimensions[0] / INCHES_PER_FOOT
		engine_center = (engine_back + engine_front) / 2

		#seats section, assuming very small seating area and a fuel tank immediately behind
		length_seats = 3 # ft
		seats_front = engine_back
		seats_back = seats_front  + length_seats
		seats_center = (seats_back + seats_front) / 2

		fuselage_width = 4.5 #ft, this seems to have been arbitrary, its scaled up from minimum for sitting
		fuselage_height = 4.5 #ft, this seems to have been arbitrary, its scaled up from minimum for sitting

		fuel_tank_front = seats_back
        # size of the fuel tank section based on the fuel volume estimate
		fuel_tank_back = fuel_tank_front + self.fuel_vol_est / ((fuselage_height - 0.5) * (fuselage_width - 0.5))
		fuel_tank_center = (fuel_tank_back + fuel_tank_front) / 2
		weight_no_wing = (self.engine.weight + self.flight.weight_stuff + self.fuel_weight_est) # stuff includes food, people, and equipment
		cog_no_wing = (engine_center * self.engine.weight + seats_center * self.flight.weight_stuff + fuel_tank_center * self.fuel_weight_est) / weight_no_wing

		wing_mac = cog_no_wing
		print((self.cr-self.ct) * self.y_bar/(self.b/2))
		print(self.y_bar)
		wing_center_geo = wing_mac + 0.25 * self.c_bar # geometric wing centrepoint
		wing_cog = wing_mac + (0.4-0.25) * self.c_bar
		wing_weight = 2.5 * self.S
		overall_cog = (cog_no_wing * weight_no_wing + wing_cog * wing_weight) / (wing_weight + weight_no_wing)
		taper_length = 4
		taper_front = fuel_tank_back
		taper_back = taper_front + taper_length
		return cog_no_wing, overall_cog, [0, engine_back, seats_back, fuel_tank_back, taper_back], wing_center_geo


# Section 8.6.5
	def calculate_tail(self):
		v_ht = 0.7 # from literature
		v_vt = 0.04 # from literature
		l_ht_from_tip = self.body_lengths[4]-1 # ft, chosen arbitrarily, but generally near total fuselage length

		S_ht = v_ht*self.c_bar*self.S/(l_ht_from_tip-self.cog_wing)

		l_vt = (l_ht_from_tip-self.cog_wing) - 1 # arbitrarily, they used 1.13, I rounded
		S_vt = v_vt * self.b*self.S / l_vt

		print(f"Tail areas: {S_ht} ft^2 (horizontal), {S_vt} ft^2 (vertical)")

		b_ht, cr_ht, ct_ht, y_bar_ht, c_bar_ht = self.determine_wing_geometry(S_ht, 4, 0.5) # arbitrary AR and taper
		b_vt, cr_vt, ct_vt, z_bar_vt_by_two, c_bar_vt = self.determine_wing_geometry(S_vt, 1.5, 0.5) # arbitrary AR and taper
		z_bar_vt = 2 * z_bar_vt_by_two

		print("\nTail dimension (span, root, tip)")
		print(f"Horizontal: {b_ht}, {cr_ht}, {ct_ht}")
		print(f"Vertical: {b_vt}, {cr_vt}, {ct_vt}")

		return S_ht, S_vt

# Sec 8.6.6
	def propeller_calc(self):
		D = 18 * (self.P/FT_LBS_PER_S_PER_HP) ** (1/4) / INCHES_PER_FOOT
		V_tip_0 = np.pi * (self.engine.rpm/60) * D
		V_tip = np.sqrt(V_tip_0**2 + self.flight.v_max**2)

		print(f"Propeller: {D} ft, {V_tip_0} ft/s at tip, {V_tip} ft/s relative")
		return D


# Sec 8.6.7
	def landing_gear_and_wing_placement_calc(self):
		height_off_ground = 0.75 + self.prop_diam/2
		Vht = 0.7 # from literature
		x_n = 0.1 * self.c_bar + self.cog_wing
		x_ac_wing = x_n - Vht
		x_geometric_center = x_ac_wing - (0.25 * self.c_bar)
		x_leading = x_geometric_center + 0.5 * self.cr
		print(f"Tried to place landing gear {x_geometric_center} ft back from nose to be on wing")
		print(f"Entire wing is forward of cog, using bicycle arrangement instead")

		A_d = 1.51
		B_d = 0.349

		A_w = 0.715
		B_w = 0.312

		nose_wheel_placement = 2 # ft, arbitrary
		main_wheel_placement = 4 # ft back from cog because we can't place it on wing, its too far forward
		x_1 = self.cog_wing - nose_wheel_placement
		x_2 = main_wheel_placement
		x_3 = x_1 + x_2

		F_M = self.weight * x_1 / x_3 # best of luck
		F_N = self.weight * x_2 / x_3
		wheel_N_diameter = A_d * F_N ** B_d
		wheel_N_width = A_w * F_N ** B_w
		wheel_M_diameter = A_d * F_M ** B_d
		wheel_M_width = A_w * F_M  ** B_w

		print(f"Wheel_N: {wheel_N_diameter} in diameter, {wheel_N_width} in width")
		print(f"Wheel_M: {wheel_M_diameter} in diameter, {wheel_M_width} in width")
		print(f"Bicycle wheel requires two extra little wheels for roll stability, these have been ignored in terms of load bearing")

# Sec 8.7
	def updated_weight_estimate(self):
		wing_weight = 1.25 * self.S
		horiz_tail_weight = 1.25 * self.Sht
		vert_tail_weight = 1.25 * self.Svt
		fraction_other = 0.05
		fraction_landing_gear = 0.057
		engine_bay_length = self.body_lengths[1] - self.body_lengths[0]
		main_fuselage_length = self.body_lengths[3] - self.body_lengths[1]
		final_taper_length = self.body_lengths[4] - self.body_lengths[3]

		Area_nose = np.pi * self.fuselage_diam * np.sqrt(engine_bay_length**2+self.fuselage_diam**2)
		Area_body = self.fuselage_diam ** 2 * main_fuselage_length
		Area_taper = np.pi * self.fuselage_diam * np.sqrt(final_taper_length**2+self.fuselage_diam**2)
		S_wet = Area_nose + Area_body + Area_taper
		fuselage_weight = S_wet * 1.4

		weight = (self.fuel_weight_est + self.flight.weight_stuff + wing_weight + horiz_tail_weight + vert_tail_weight + fuselage_weight + self.engine.weight) / (1-fraction_other-fraction_landing_gear)
		print(f"Updated weight estimate: {weight} lbs")
		print(f"Updated weight changed by: {weight - self.weight}")
		return weight

	def other_configuration_stuff(self):
		print("5 degree dihedral from previous designs")
		print("make the control surfaces 30% of local chord")
		print("hatch located next to seats")
		print("make the airlerons 50% of the wing length")


    # calculating how much fuel *could* be stored in the wings, to assess viability
	def fuel_in_wings(self):
		geo = self.airfoil.airfoil_geometry
		N_span = 100 # number of slices in the spanwise direction (for single wing)
		d_span = (self.b/2) / N_span
		N_chord = len(self.airfoil.airfoil_geometry)

		volume = 0
		length = 0
		for span_index in range(N_span):
			c = self.cr - (self.cr-self.ct) * ((span_index+0.5)/N_span)
			length = length+d_span
			for chord_index in range(N_chord-1):
				slice_width = (geo[chord_index+1,0] - geo[chord_index,0]) * c
				average_top = (geo[chord_index+1,1]+geo[chord_index,1])/2 * c
				average_bottom = (geo[chord_index+1,2]+geo[chord_index,2])/2 * c
				slice_height = average_top - average_bottom
				volume = volume + slice_height * slice_width * d_span
		return volume

# TODO: relocate Flight and Engine into the above sections where they make the most sense, I just don't really know where they go right now
class Flight:
	def __init__(self, range, crew, time, safety, stall):
		self.range = range * KM_TO_MILES * MILES_TO_FEET
		self.crew = crew
		self.time = time
		self.weight_stuff = self.crew * (CREW_MEMBER_WEIGHT + CREW_PAYLOAD_WEIGHT) * KG_TO_LBS
		self.safety = safety
		self.v_max = 200 * KM_TO_MILES * MILES_TO_FEET / SECONDS_PER_HOUR
		self.v_stall = stall
		self.v_flare = self.v_stall * 1.23
		self.landing_distance = 5280 * 3 # ft, assuming a ~3 mile runway is about the max practical (assuming access to a military base)
		self.takeoff_distance = 5280 * 3 # ft
		self.obstacle_height = 10 # ft
		self.rc_max = 2 #ft/s

class Engine:
	def __init__(self, filename):
		with open(filename) as f:
			data = json.load(f)
		self.nominal_power = data['power'] # we are going to use nominal power, because we can just fly at 18000 feet so turbocharging doesn't matter
		self.turbocharged_to = data['turbocharged_to']
		self.weight = 1.4 * data['nominal_weight'] # 1.4 is estimate from raynor for installed weight
		self.dimensions = [data['length'], data['width'], data['height']]
		self.rpm = data['rpm']

def main():
	# Flight information
	R = 40000 # km
	CREW = 2 # people
	FLIGHT_TIME = 10 # days
	SAFETY_FACTOR = 1.01
	V_STALL = 180 * KM_TO_MILES * MILES_TO_FEET / SECONDS_PER_HOUR # this is already pretty high, but we can mess with it

	# this just stores the goal, vs details of aircraft
	flight = Flight(R, CREW, FLIGHT_TIME, SAFETY_FACTOR, V_STALL)

	# could be swapped for a different airfoil profile if we wanted
	airfoil = Airfoil('./xf-r1046-il-1000000.csv', './airfoil_profile.csv')

	# engine details
	engine = Engine('./textron_540_v.json')


	# creating the aircraft sets up an initial guess for key quantities
	aircraft = Aircraft(flight, airfoil, engine)

	airfoil.plot_airfoil()

if __name__ == '__main__':
	main()