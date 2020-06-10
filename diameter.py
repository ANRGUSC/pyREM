import numpy as np 
import math
import scipy
from scipy.optimize import root

D_0 = 1.00*10**-5 #initial diameter of droplet
RELATIVE_HUMMIDITY = 60 #relative hummidity
TEMPERATURE = 293.15 # ambient temperature in Kelvin
RV = 461.52 #J/kgK specific gas constant for water
RHO_D = 1000 #density of droplet in kg/m^3; same as RHO_P
RHO_A = 1.21 #also called RHO, density of air in kg/m^3 
RHO_P = RHO_D
VISCOSITY = 1.81*10**-5 #viscosity of air in Pa s


def diameter_polynomial(time,initial_D=D_0):
    molec_diff = (2.16*10**-5)*(TEMPERATURE/273.15)**1.8 #molecular diffusivity of water vapor
    p_sat = 611.21*math.exp((19.843-(TEMPERATURE/234.5))*((TEMPERATURE-273.15)/(TEMPERATURE-16.01)))
    p_infin = p_sat*RELATIVE_HUMMIDITY/100

    k = ((8*molec_diff*(p_sat-p_infin)*(initial_D**2)*time)/(RHO_P*RV*TEMPERATURE))
    m = -initial_D**2
    p = np.poly1d([1, 0, m, 0, k])

    roots = max(np.roots(p))
    d = roots

    if np.iscomplex(d) == True:
       d = 7.00*10**-6

    return d

 
def velocity_exponential(time,initial_D=D_0):
    if time <= 0:
        return (RHO_P*initial_D**2*G)/(18*math.pi*VISCOSITY) #Stoke's Law for small velocities

    d = diameter_polynomial(time,initial_D) 
    n = 10.8*VISCOSITY*((RHO_A*d)/VISCOSITY)**0.687 
    p = 4*(d**2)*(RHO_D-RHO_A) 
    m = 72*VISCOSITY

    roots = root(lambda v: n*v**(2.687)+m*v**2-p*v,0.1)
    
    return roots.x[0]

if __name__ == '__main__':
    #diameter_polynomial(0.01)
    velocity_exponential(60)
