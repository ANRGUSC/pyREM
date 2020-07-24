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
D_0 = 1.00*10**-5 #initial diameter of droplet in micro meters
A = 0.06 #given constant in dispersion coefficient equation
B = 0.92 #given constant in dispersion coefficient equation
NUMBER_OF_DROPLETS = 1 #number of droplets emitted (q)
X_0 = 0 #initial horizontal position
Z_0 = 0 #initial vertical position
RESPIRATORY_RATE = 0.25 #avg number of breaths/second 
V_X = 1 #horizontal velocity of air surrounding the droplets in m/s
RELATIVE_HUMMIDITY = 60 #default relative hummidity
TEMPERATURE = 293.15 #default ambient temperature in Kelvin
X_AWAY = 2 #default distance 2 meters away from source 


def diameter_polynomial(time,r_h,initial_D):
    '''This function estimates the droplet's diameter in micrometers 
    as a function of time and initital size by finding the real roots of 
    the diameter polynomial. If the roots are complex, the droplet diameter has 
    reached its minimum, dmin, and is estimated as a percentage of the initial diameter. 

    Parameters:
        time (float): time at which the droplet diameter will be calculated
        r_h (int): relative humidity 
        initial_D (float): initial droplet size in micrometers
                    
    Returns:
        d (float): Returns d, a float value representing the diameter 
        of the droplet after t seconds. 
    '''
    molec_diff = (2.16*10**-5)*(TEMPERATURE/273.15)**1.8 #molecular diffusivity of water vapor
    p_sat = 611.21*math.exp((19.843-(TEMPERATURE/234.5))*((TEMPERATURE-273.15)/(TEMPERATURE-16.01)))
    p_infin = p_sat*r_h/100

    k = ((8*molec_diff*(p_sat-p_infin)*(initial_D**2)*time)/(RHO_P*RV*TEMPERATURE))
    m = -initial_D**2
    p = np.poly1d([1, 0, m, 0, k])

    roots = max(np.roots(p))
    d = roots

    if np.iscomplex(d) == True:
        d = 0.71*initial_D

    return d

def terminal_velocity(time,r_h,initial_D):
    ''' This function estimates the terminal velocity in m/s of the droplet as a function of time 
    and initial droplet size. For small velocities, v_t is calculated using Stoke's Law. Otherwise,
    it is calculated by finding the roots of the velocity exponential.

    Parameters:
        time (float): time at which the terminal velocity will be calculated
        r_h (int): relative humidity 
        initial_D (float): initial droplet size in micrometers

    Returns: 
        v_t (float): v_t, a float value representing the terminal velocity of the droplet 
        after t seconds.
    '''
    if time <= 0:
        v_t = (RHO_P*initial_D**2*G)/(18*math.pi*VISCOSITY) #Stoke's Law for small velocities
    else:
        d = diameter_polynomial(time,r_h,initial_D) 
        n = 10.8*VISCOSITY*((RHO_A*d)/VISCOSITY)**0.687 
        p = 4*(d**2)*(RHO_D-RHO_A)*G 
        m = 72*VISCOSITY
        roots = root(lambda v: n*v**(2.687)+m*v**2-p*v,1)
        v_t = roots.x[0]

    return v_t


def position(time,r_h,initial_D): 
    ''' This function estimates the horizontal and vertical distance of droplet after t seconds. 
    The vertical distance, z_d, is calculated using an integral since the terminal velocity continues 
    to change until the droplet's diameter reaches its minimum, dmin.

    Parameters:
        time (float): time at which the x_d and z_d values are calculated
        r_h (int): relative humidity
        initial_D (float): initial droplet size in micrometers

    Returns: 
        (x_d,z_d): a 2-tuple of float values containing the x and z positions of the droplet in meters
        after t seconds.
    ''' 
    if time <= 0:
        return (X_0, Z_0)

    d = diameter_polynomial(time,r_h,initial_D)
    v_t = terminal_velocity(time,r_h,initial_D)

    v_integral = integrate.quad(terminal_velocity, 0, time, args=(r_h,initial_D,))

    x_d = X_0 + V_X*time
    z_position = Z_0-v_integral[0]

    if z_position >= -2:
        z_d = z_position 
    else:
        z_d = -2 #droplet reaches the ground

    distance_tuple = (x_d,z_d)

    return distance_tuple

def concentration(time,r_h,initial_D): 
    ''' Each breath is modeled as an ever expanding Gaussian distribution or puff containing
    thousands of respiratory droplets. This function estimates the concentration of the puff at a 
    particular time. 

    Parameters:
        time (float): time in seconds
        r_h (int): relative humidity
        initial_D (float): initial droplet size in micrometers

    Returns: 
        conc_of_puff (float): a float value representing the concentration of the puff that interacts with a 
        person 2 meters from the infected source.
    ''' 
    distance_tuple = position(time,r_h,initial_D)
    x_d = distance_tuple[0]
    z_d = distance_tuple[1]
    sigma = A*(x_d**B) #dispersion coefficient 
    conc_of_puff = (NUMBER_OF_DROPLETS/((math.sqrt(2*math.pi)*sigma))**3)*math.exp((-1/(2*sigma**2))*((X_AWAY-x_d)**2+z_d**2))

    return conc_of_puff

def exposure_per_breath(time,r_h,initial_D): 
    ''' This function estimates the dose of respiratory droplets that a person is exposed to by 
    integrating the puff over time. The function uses a numerical approach to calculate the integral
    by using the trapezoidal rule.  

    Parameters:
        time (float): time in seconds
        r_h (int): relative humidity
        initial_D (float): initial droplet size in micrometers

    Returns: 
        exposure (float): A float value representing the exposure per breath 
        or concentration of the puff integrated over time.
    ''' 
    time_array = []
    conc_array = []
    for i in range(1,501,1):
      tm = (i)/(100.0)
      conc = concentration(tm,r_h,initial_D)
      conc_array.append(conc)
      time_array.append(tm)
    
    y = conc_array
    x = time_array
    for i in range(0,len(x)):
       x[i] = round(x[i],5)
       y[i] = round(y[i],5)
    
    exposure = np.trapz(y,x)

    return exposure 

    
def total_exposure(time,r_h=RELATIVE_HUMMIDITY,initial_D=D_0):
    ''' This function estimates the total exposure by multiplying the 
    exposure per breath by the number of breaths taken in t seconds.

    Parameters:
        time (float): time in seconds
        r_h (int): relative humidity
        initial_D (float): initial droplet size in micrometers

    Returns: 
        total_dosage (float): a float value representing the total dosage a person 
        is exposed to after several breaths are taken from an infected source. 
    '''
    exposure_p_breath = exposure_per_breath(time,r_h,initial_D)
    number_of_breaths = RESPIRATORY_RATE*time
    total_dosage = exposure_p_breath*number_of_breaths

#    print(total_dosage)

    return total_dosage
    
#example usage, for testing
if __name__ == '__main__':
#    total_exposure(5)

'''
#   plot for concentration vs droplet size with humidity varied
    t = 5
    initial_D_list = list(np.arange(1*10**-6, 1*10**-4, 2*10**-6))
    hummidity = [50,60,70,80]
    for r in hummidity:
        exposure_array = []
#        print(r)
        for init_D in initial_D_list:
#            print(init_D)
            exposure = exposure_per_breath(t,r,init_D)/(12.3487382)
            exposure_array.append(exposure)
        plt.plot(initial_D_list,exposure_array, label = "RH = " + str(r))
    plt.xlabel('Droplet Size')
    plt.ylabel('Concentration of Droplets')
    plt.title('Concentration vs Droplet Size Graph')
    plt.legend()
    plt.show() 
'''
