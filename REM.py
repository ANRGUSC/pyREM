import math
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root

RHO_A = 1.21 #density of air in kg/m^3 
RHO_D = 1000 #density of droplet in kg/m^3
RHO = RHO_A
RHO_P = RHO_D
G = 9.81 #gravitational acceleration in m/s^2
VISCOSITY = 1.81*10**-5 #viscosity of air in Pa s
RV = 461.52 #J/kgK specific gas constant for water
D_0 = 1.0*10**-5 #initial diameter of droplet in micro meters
A = 0.06 #given constant in dispersion coefficient equation
B = 0.92 #given constant in dispersion coefficient equation
NUMBER_OF_DROPLETS = 1 #number of droplets emitted (q)
X_0 = 0 #initial horizontal position
Z_0 = 0 #initial vertical position
RESPIRATORY_RATE = 0.25 #avg number of breaths/second 
V_X = 1 #horizontal velocity of the air surrounding the droplets in m/s
RELATIVE_HUMIDITY = 60 #default relative humidity
TEMPERATURE = 293.15 #default ambient temperature in Kelvin
X_AWAY = 2 #default distance 2 meters away from source 


def diameter_polynomial(time,temp,r_h,initial_D): 
    '''This function estimates the droplet's diameter in micrometers  
    by finding the real roots of the diameter polynomial. If the roots are complex, 
    the droplet diameter has reached its minimum, dmin, and is estimated at time = t_crit,
    where the discrimiant of the polynomial is zero.

    Parameters:
        time (float): time at which the droplet diameter will be calculated
        temp (float): ambient temperature in Kelvin 
        r_h (int): relative humidity 
        initial_D (float): initial droplet size in micrometers
                    
    Returns:
        d (float): Returns d, a float value representing the diameter 
        of the droplet after t seconds. 
    '''
    molec_diff = (2.16*10**-5)*(temp/273.15)**1.8 # molecular diffusivity of water vapor
    p_sat = 611.21*math.exp((19.843-(temp/234.5))*((temp-273.15)/(temp-16.01))) # saturation water vapor pressure
    p_infin = p_sat*r_h/100 # ambient water vapor pressure
    t_crit = (RHO_P*RV*temp*(initial_D**2))/(32*molec_diff*(p_sat-p_infin)) # time when Discriminant is 0

    k = ((8*molec_diff*(p_sat-p_infin)*(initial_D**2)*time)/(RHO_P*RV*temp))
    m = -initial_D**2
    p = np.poly1d([1, 0, m, 0, k])
    roots = max(np.roots(p))
    
    if time <= t_crit:
        d = roots
    else:
       d = diameter_polynomial(t_crit,temp,r_h,initial_D)
#    d = roots

#    if np.iscomplex(d) == True:
#        d = 0.71*initial_D

    return d


def terminal_velocity(time,temp,r_h,initial_D):
    '''This function estimates the terminal velocity in m/s of the droplet as a function of time, 
    temperature, humidity and initial droplet size. For small velocities, v_t is calculated 
    using Stoke's Law. Otherwise, it is calculated by finding the roots of the velocity exponential.

    Parameters:
        time (float): time at which the terminal velocity will be calculated
        temp (float): ambient temperature in Kelvin 
        r_h (int): relative humidity 
        initial_D (float): initial droplet size in micrometers

    Returns: 
        v_t (float): v_t, a float value representing the terminal velocity of the droplet 
        after t seconds
    '''
    if time <= 0:
        v_t = (RHO_P*initial_D**2*G)/(18*math.pi*VISCOSITY) #Stoke's Law for small velocities
    else:
        d = diameter_polynomial(time,temp,r_h,initial_D) 
        n = 10.8*VISCOSITY*((RHO_A*d)/VISCOSITY)**0.687 
        p = 4*(d**2)*(RHO_D-RHO_A)*G 
        m = 72*VISCOSITY
        roots = root(lambda v: n*v**(2.687)+m*v**2-p*v,1)
        v_t = roots.x[0]

    return v_t

def position(time,temp,r_h,initial_D): 
    ''' This function estimates the horizontal and vertical position of droplet after t seconds. 
    The vertical distance, z_d, is calculated using an integral since the terminal velocity continues 
    to change until the droplet's diameter reaches its minimum, dmin.

    Parameters:
        time (float): time at which the x_d and z_d values are calculated
        temp (float): ambient temperature in Kelvin 
        r_h (int): relative humidity 
        initial_D (float): initial droplet size in micrometers

    Returns: 
        (x_d,z_d): a 2-tuple of float values containing the x and z positions of the droplet in meters
    '''	
    if time <= 0:
        return (X_0, Z_0)

    v_t = terminal_velocity(time,temp,r_h,initial_D)
    v_integral = integrate.quad(terminal_velocity, 0, time, args=(temp,r_h,initial_D,))

    x_d = X_0 + V_X*time 
    z_position = Z_0-v_integral[0]

    if z_position >= -2:
        z_d = z_position 
    else:
        z_d = -2 #droplet reaches the ground

    distance_tuple = (x_d,z_d)

    return distance_tuple

def concentration(time,x_away,temp,r_h,initial_D): 
    ''' Each breath is modeled as an expanding Gaussian puff containing
    thousands of respiratory droplets. This function estimates the concentration of the puff at a 
    particular time. 

    Parameters:
        time (float): time in seconds 
        x_away (float): distance x meters away from an infected source 
        temp (float): ambient temperature in Kelvin 
        r_h (int): relative humidity 
        initial_D (float): initial droplet size in micrometers

    Returns: 
        conc_of_puff (float): a float value representing the concentration of the puff that interacts with a 
        person x meters from an infected source.
    '''	
    distance_tuple = position(time,temp,r_h,initial_D)
    x_d = distance_tuple[0]
    z_d = distance_tuple[1]
    sigma = A*(x_d**B) #dispersion coefficient 
    conc_of_puff = (NUMBER_OF_DROPLETS/((math.sqrt(2*math.pi)*sigma))**3)*math.exp((-1/(2*sigma**2))*((x_away-x_d)**2+z_d**2))

    return conc_of_puff

def exposure_per_breath(time,x_away,temp,r_h,initial_D): 
    '''This function estimates the dose of respiratory droplets that a person is exposed to by 
    integrating the puff over time. The function uses the quad function to calculate the integral
    using 50 subdivisions.  

    Parameters:
        time (float): time in seconds that represents the upper limit of the integral
        x_away (float): distance x meters away from the infected source 
        temp (float): ambient temperature in Kelvin 
        r_h (int): relative humidity 
        initial_D (float): initial droplet size in micrometers

    Returns: 
        exposure (2-tuple float): A 2-tuple of float value containing the 
        concentration of the puff integrated over time and the possible 
        numerical error in the integrand from the use of quad 
    ''' 
    exposure = integrate.quad(concentration, 0, time, args=(x_away,temp,r_h,initial_D,), limit=50) #integrating with respect to time

    return exposure

def total_exposure(time,x_away=X_AWAY,temp=TEMPERATURE,r_h=RELATIVE_HUMIDITY, initial_D=D_0):
    '''This function estimates the total exposure by multiplying the 
    exposure per breath by the number of breaths taken in t seconds.

    Parameters:
        time (float): time in seconds
        x_away (float): proximity set to the default value of 2 meters
        temp (float): temperature set to the default value of 293.15 K (20 C)
        r_h (int): humidity set to the default value of 60
        initial_D (float): initial droplet size set to the default value of 10 um

    Returns: 
        total_dosage (float): a float value representing the total dosage a person 
        is exposed to after several breaths are taken from an infected source. 
    '''
    exposure_tuple = exposure_per_breath(time,x_away,temp,r_h,initial_D)
    number_of_breaths = RESPIRATORY_RATE*time
    total_dosage = exposure_tuple[0]*number_of_breaths

#    print(total_dosage)

    return total_dosage


#example usage, for testing
if __name__ == '__main__':
    total_exposure(5) #total accumulated exposure after 5 seconds 
    
