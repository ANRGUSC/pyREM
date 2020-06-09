import numpy as np 
import math
import scipy
from scipy import optimize
from scipy.optimize import fsolve

D_0 = 1.00*10**-5 #initial diameter of droplet
RELATIVE_HUMMIDITY = 60 #relative hummidity
TEMPERATURE = 293.15 # ambient temperature in Kelvin
RV = 462.52 #J/kgK specific gas constant for water
RHO_D = 1000 #density of droplet in kg/m^3; same as RHO_P
RHO_A = 1.21 #also called RHO, density of air in kg/m^3 
RHO_P = RHO_D
VISCOSITY = 1.81*10**-5 #viscosity of air in Pa s
#N = 10.8*VISCOSITY*((RHO_A*D_0)/VISCOSITY)**0.687 #use d  from function not D_0
#P = 4*(D_0**2)*(RHO_D-RHO_A) #use d from function not D_0
#M = 72*VISCOSITY

def diameter_polynomial(time,initial_D=D_0): #needs to return d 
    molec_diff = (2.16*10**-5)*(TEMPERATURE/273.15)**1.8 #molecular diffusivity of water vapor
    p_sat = 611.21*math.exp((19.843-(TEMPERATURE/234.5))*((TEMPERATURE-273.15)/(TEMPERATURE-16.01)))
    p_infin = p_sat*RELATIVE_HUMMIDITY/100

    k = ((8*molec_diff*(p_sat-p_infin)*(initial_D**2)*time)/(RHO_P*RV*TEMPERATURE))
    m = initial_D**2

    p = np.poly1d([1, 0, m, 0, k])
    print(np.roots(p))
    #print(k)
    #print(m)
    return 

def droplet_diameter(time, initial_D): 
    ''' This function estimates the droplet's diameter in micrometers as a function of time (s) that also depends
     on the relative hummidity (RH) and the temperature (T) in Kelvin which are set as constants in this program.

    Parameters:
        time (float): This parameter represents the time at which the diameter of the droplet will be calclulated, 
                    since the size of the droplet evaporates over time.
    Returns:
        d (float): Returns d, a float value representing the diameter of the droplet after it has 
                evaporated after t seconds. 
    '''
    if time <= 0:
        return initial_D

    molec_diff = (2.16*10**-5)*(TEMPERATURE/273.15)**1.8 #molecular diffusivity of water vapor
    p_sat = 611.21*math.exp(((19.843-TEMPERATURE)/234.5)*((TEMPERATURE-273.15)/(TEMPERATURE-16.01)))
    p_infin = p_sat*RELATIVE_HUMMIDITY/100
    beta = (8*molec_diff*(p_sat-p_infin))/((initial_D**2)*RHO_P*RV*TEMPERATURE) #evaporation rate
    #beta = (8*molec_diff*(p_sat-p_infin))/((droplet_diameter(time-0.005, D_0)**2)*RHO_P*RV*TEMPERATURE) #using previous droplet size

    d_min = 0.44*initial_D
    d = max(initial_D*math.sqrt(max(1-beta*time,0)), d_min)

    return d

D = droplet_diameter(0.01,D_0)
N = 10.8*VISCOSITY*((RHO_A*D)/VISCOSITY)**0.687 
P = 4*(D**2)*(RHO_D-RHO_A) 
M = 72*VISCOSITY

def v(x):
	return N*x**2.687+M*x**2-P*x
 
def velocity_polynomial(time,initial_D=D_0):
    d = droplet_diameter(time,initial_D) #will replace this with diameter_polynomial after it is working 
    n = 10.8*VISCOSITY*((RHO_A*d)/VISCOSITY)**0.687 
    p = 4*(d**2)*(RHO_D-RHO_A) 
    m = 72*VISCOSITY
    #function = n*x**2.687+m*x**2-p*x
    return 

if __name__ == '__main__':
    #diameter_polynomial(0.01)
    roots = fsolve(v,0)
    print(roots)







