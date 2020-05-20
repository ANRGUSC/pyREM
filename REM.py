import time
import math
import sympy as sp


RHO_A = 1.21 #density of air in kg/m^3
RHO_D = 1000 #density of droplet in kg/m^3
G = 9.81 #gravitational acceleration in m/s^2
VISCOSITY = 1.81*10**-5 #viscosity of air in Pa s
RV = 462.52 #J/kgK specific gas constant for water
D_0 = 1.00*10**-5 #initial diameter of droplet

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
    v_x = math.sqrt((4*d*(RHO_D - RHO_A)*G)/3*RHO_A*drag_coef)
    print(v_x)
    return;



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

    print(d)
    return;

def x_position(v_x,v_t):
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

if __name__ == '__main__':
    #terminal_velocity(0.00005,1)
    #droplet_diameter(0.04,60,230)



