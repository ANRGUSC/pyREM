import time
import math
import scipy.integrate as integrate
import matplotlib.pyplot as plt

RHO_A = 1.21 #density of air in kg/m^3
RHO_D = 1000 #density of droplet in kg/m^3
G = 9.81 #gravitational acceleration in m/s^2
VISCOSITY = 1.81*10**-5 #viscosity of air in Pa s
RV = 462.52 #J/kgK specific gas constant for water
D_0 = 1.00*10**-5 #initial diameter of droplet
A = 0.06 #given constant in dispersion coefficient equation
B = 0.92 #given constant in dispersion coefficient equation
NUMBER_OF_DROPLETS = 1 #number of droplets emitted (q)
X_0 = 0 #initial X position
Z_0 = 2 #Initial vertical position
RESPIRATORY_RATE = 0.25 #breaths/second from avg of 15 bpm
RELATIVE_HUMMIDITY = 60 #relative hummidity
TEMPERATURE = 293.15 # ambient temperature in Kelvin
V_X = 1 #horizontal velocity 1m/s
X_AWAY = 4 #a distance X meters away from source 

def terminal_velocity(d):
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

    reynolds_p = RHO_A*V_X*d/VISCOSITY #reynolds number calculation
    drag_coef = 24*(1+0.15*reynolds_p**0.687)/reynolds_p #Drag Coefficient
    v_t = math.sqrt((4*d*(RHO_D - RHO_A)*G)/3*RHO_A*drag_coef)
    return v_t


def droplet_diameter(time): 
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
    molec_diff = (2.16*10**-5)*(TEMPERATURE/273.15)**1.8 #molecular diffusivity of water vapor
    p_sat = 611.21**(((19.843-TEMPERATURE)/234.5)*((TEMPERATURE-273.15)/(TEMPERATURE-16.01)))
    p_infin = p_sat*RELATIVE_HUMMIDITY/100
    beta = (8*molec_diff*(p_sat-p_infin))/((D_0**2)*RV*TEMPERATURE) #evaporation rate

    d_min = 0.44*D_0
    d = max(D_0*math.sqrt(max(1-beta*time,0)), d_min)

    return d

def position(time): # function that returns x and z position as tuple
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
    d = droplet_diameter(time)
    v_t = terminal_velocity(d)
    time_to_hit_ground = Z_0/v_t
    x_d = X_0 + V_X*min(time, time_to_hit_ground)
    z_d = max(Z_0-v_t*time,0) #take into account droplet hitting the ground

    distance_tuple = (x_d,z_d)
    return distance_tuple

def concentration(time,x_away=X_AWAY): 
    distance_tuple = position(time)
    sigma = A*(x_away**B)
    x_d = distance_tuple[0]
    z_d = distance_tuple[1]
    conc_of_puff = (NUMBER_OF_DROPLETS/(math.sqrt(2*math.pi*sigma))**3)*math.exp(((-1/2*sigma**2)*((x_away-x_d)**2)+z_d**2))
    return conc_of_puff

def exposure_per_breath(time,x_away=X_AWAY): 
    exposure = integrate.quad(concentration, 0, time, args=(x_away,)) #keep x_away constant while integrating
    #print(exposure)
    return exposure

def total_exposure(time=5,x_away=X_AWAY):
    exposure_tuple = exposure_per_breath(time,x_away)
    number_of_breaths = RESPIRATORY_RATE*time
    total_dosage = exposure_tuple[0]*number_of_breaths
    #print(total_dosage)
    return total_dosage
    


if __name__ == '__main__':
    #exposure_per_breath(5,4)
    #total_exposure(5)
    #exposure_over_distance(50,4)
    exposure_array = []
    for i in range(0,11):
        exposure = total_exposure(i)
        exposure_array.append(exposure)
    #print(exposure_array)
    droplet_size_array = [0,10,20,30,40,50,60,70,80,90,100]
    plt.plot(droplet_size_array,exposure_array)
    plt.xlabel('Droplet Size')
    plt.ylabel('Concentration of Droplets')
    plt.title('Concentration vs Droplet Size Graph')
    plt.show()
    #multiple plots example
    #plt.plot(diameter_array,exposure_array1,
           #diameter_array,exposure_array2)
