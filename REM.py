import time
import math
from scipy.integrate import quad


RHO_A = 1.21 #density of air in kg/m^3
RHO_D = 1000 #density of droplet in kg/m^3
G = 9.81 #gravitational acceleration in m/s^2
VISCOSITY = 1.81*10**-5 #viscosity of air in Pa s
RV = 462.52 #J/kgK specific gas constant for water
D_0 = 1.00*10**-5 #initial diameter of droplet
A = 0.06 #given constant in dispersion coefficient equation
B = 0.92 #given constant in dispersion coefficient equation
#NUMBER_OF_DROPLETS =  #number of droplets emitted (q)

def terminal_velocity(d,v_x):
    ''' This function estimates the terminal velocity (units um/s ?) of a respitory droplet as a function of
    the droplet's diamter "d" and horizontal velocity "v_x"

    Parameters:
        d (float): A float representing the the diamter of the droplet which will be used to calculate the 
                Reynolds number and then used to find the drag coefficient "CD"
        v_x (float): In this model, the horizontal velocity of the droplet is equivalent to the velocity of the 
                of the surrounding air due to the small size of the droplet. In most cases, v_x is assumed 
                to be 1 m/s.

    Returns: 
        v_t (float): This function returns the terminal velocity of the droplet, v_t, when Fdrag = Fgrav

    '''

    reynolds_p = RHO_A*v_x*d/VISCOSITY #reynolds number calculation
    drag_coef = 24*(1+0.15*reynolds_p**0.687)/reynolds_p #Drag Coefficient
    v_t = math.sqrt((4*d*(RHO_D - RHO_A)*G)/3*RHO_A*drag_coef)
    #print(v_t)
    return v_t


def droplet_diameter(time,r_h,temp):
    ''' This function estimates the droplet's diameter in meters, as a function of time (s), the relative hummidity (RH),
    and the temperature, T in Kelvin.

    Parameters:
        time (float): This parameter represents what time the dimater of the droplet will be calclulated, since
                    the size of the droplet evaporates over time.
        r_h (int): An integer repsenting the relative hummidity of the ambient conditions. 
        temp (float): This parameter represents the temperature of the ambient conditions in Kelvin which will be used
                     in the evaporation rate, beta to solve for the diameter. 

    Returns:
        d (float): This function returns d, a value representing the diameter of the droplet after it has 
                evaporated after t seconds. 

    '''
    molec_diff = (2.16*10**-5)*(temp/273.15)**1.8 #molecular diffusivity of water vapor
    p_sat = 611.21**(((19.843-temp)/234.5)*((temp-273.15)/(temp-16.01)))
    p_infin = p_sat*r_h/100
    beta = (8*molec_diff*(p_sat-p_infin))/((D_0**2)*RV*temp) #evaporation rate

    d_min = 0.44*D_0
    d = max(D_0*math.sqrt(1-beta*time), d_min)

    #print(d)
    return d

def x_position(x_0,z_0,v_x,v_t,time):
    ''' This function estimates the horizontal distance the droplet has travelled at a certain terminal velocity and
    horizontal velocity.  

    Parameters:
        v_x (float): In this model, the horizontal velocity of the droplet is equivalent to the velocity of the 
                of the surrounding air due to the small size of the droplet. In most cases, v_x is assumed 
                to be 1 m/s.
        v_t (float): A float that corresponds to the terminal velocity of the droplet

    Returns: 
        x_d (float): returns the horizontal distance of the drop in meters.  
    '''	
    time_to_hit_ground = z_0 / v_t
    return x_0 + v_x*min(time, time_to_hit_ground)

def concentation(x_away,v_t,time):
    sigma = A*(x_away**B)
    x_d = 0  #assuming source is at position 0
    z_d = v_t*time - 0.5*G*time**2
    conc_of_puff = (NUMBER_OF_DROPLETS/(math.sqrt(2*pi*sigma))**3)**((-1/2*sigma**2)*((x_away-x_d)**2)+z_d**2)
    #print(conc_of_puff)
    return conc_of_puff

def integrand(x_away,time):
	return (NUMBER_OF_DROPLETS/(math.sqrt(2*pi*sigma))**3)**((-1/2*sigma**2)*(x_away**2)+z_d**2)

def exposure_per_breath(time): 
    exposure, err = quad(integrand, 0, time)
    #print(exposure)
    return exposure

def total_exposure(num_breaths,time):
	exposure = exposure_per_breath(time) #exposure per breath
    total_dosage = exposure*num_breaths
    print(total_dosage)
    return;

if __name__ == '__main__':
    #terminal_velocity(0.00005,1)
    #droplet_diameter(0.04,60,230)
    #x_position(0,1,0.005)
    exposure_per_breath(5) #avearge breath is taken every 5s
    total_exposure(2,5)
