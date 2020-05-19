import time 
import math
import sympy as sp


CONST_rho_a = 1.21 #density of air in kg/m^3
CONST_rho_d = 1000 #density of droplet in kg/m^3
CONST_g = 9.81 #gravitational acceleration in m/s^2
CONST_viscosity = 1.81*10**-5 #viscosity of air in Pa s


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
	RE_p = CONST_rho_a*v_x*d/CONST_viscosity #reynolds number calculation
	CD = 24*(1+0.15*RE_p**0.687)/RE_p
	v_x = math.sqrt((4*d*(CONST_rho_d - CONST_rho_a)*CONST_g)/3*CONST_rho_a*CD)
	print(v_x)
	return;

if __name__ == '__main__':
	terminal_velocity(0.00005,1)
