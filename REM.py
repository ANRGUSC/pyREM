import time 
import math
import sympy as sp


rho_a = 1.21 #density of air in kg/m^3
rho_d = 1000 #density of droplet in kg/m^3
g = 9.81 #gravitational acceleration in m/s^2
viscosity = 1.81*10**-5 #viscosity of air in Pa s

def terminal_velocity(d,v_x):
	''' This function estimates the terminal velocity (units um/s ?) of a respitory droplet as a function of
	the droplet's diamter "d" and horizontal velocity "v_x"

	Parameters:
		d (float): An integer representing the the diamter of the droplet which will be used to calculate the 
				Reynolds number and then used to find the drag coefficient "CD"
		v_x (float): In this model, the horizontal velocity of the droplet is equivalent to the velocity of the 
				of the surrounding air due to the small size of the droplet. In most cases, v_x is assumed 
				to be 1 m/s.

	Returns: 
		v_t (float): An integer that returns the terminal velocity of the droplet when the drag force = Fgrav

	'''
	RE_p = rho_a*v_x*d/viscosity #reynolds number calculation
	CD = 24*(1+0.15*RE_p**0.687)/RE_p #Drag Coefficient
	v_x = math.sqrt((4*d*(rho_d - rho_a)*g)/3*rho_a*CD)
	print(v_x)
	return;

Rv = 462.52 #J/kgK specific gas constant for water
d_0 = 1.00*10**-5 #initial diameter of droplet

def droplet_diameter(time,RH,T):
	''' This function estimates the terminal velocity (units um/s ?) of a respitory droplet as a function of
	the droplet's diamter "d" and horizontal velocity "v_x"

	Parameters:
		d (float): An integer representing the the diamter of the droplet which will be used to calculate the 
				Reynolds number and then used to find the drag coefficient "CD"
		v_x (float): In this model, the horizontal velocity of the droplet is equivalent to the velocity of the 
				of the surrounding air due to the small size of the droplet. In most cases, v_x is assumed 
				to be 1 m/s.

	Returns: 
		v_t (float): An integer that returns the terminal velocity of the droplet when the drag force = Fgrav

	'''
	D = (2.16*10**-5)*(T/273.15)**1.8 #molecular diffusivity of water vapor
	P_sat = 611.21**(((19.843-T)/234.5)*((T-273.15)/(T-16.01)))
	P_infin = P_sat*RH/100
	beta = (8*D*(P_sat-P_infin))/((d_0**2)*Rv*T) #evaporation rate

	d_min = 0.44*d_0
	d = max(d_0*math.sqrt(1-beta*time), d_min)
	
	print(d)
	return;

def position(v_x,v_t):
	''' This function estimates the terminal velocity (units um/s ?) of a respitory droplet as a function of
	the droplet's diamter "d" and horizontal velocity "v_x"

	Parameters:
		d (float): An integer representing the the diamter of the droplet which will be used to calculate the 
				Reynolds number and then used to find the drag coefficient "CD"
		v_x (float): In this model, the horizontal velocity of the droplet is equivalent to the velocity of the 
				of the surrounding air due to the small size of the droplet. In most cases, v_x is assumed 
				to be 1 m/s.

	Returns: 
		v_t (float): An integer that returns the terminal velocity of the droplet when the drag force = Fgrav

	'''	

if __name__ == '__main__':
	#terminal_velocity(0.00005,1)
	droplet_diameter(0.04,60,230)



