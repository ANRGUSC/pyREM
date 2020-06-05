import math
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np

RHO_A = 1.21 #density of air in kg/m^3 
RHO_D = 1000 #density of droplet in kg/m^3; same as RHO_P
RHO = RHO_A
RHO_P = RHO_D

G = 9.81 #gravitational acceleration in m/s^2
VISCOSITY = 1.81*10**-5 #viscosity of air in Pa s
RV = 461.52 #J/kgK specific gas constant for water
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
    p_sat = 611.21*math.exp((19.843-TEMPERATURE/234.5)*((TEMPERATURE-273.15)/(TEMPERATURE-16.01)))
    p_infin = p_sat*RELATIVE_HUMMIDITY/100
    beta = (8*molec_diff*(p_sat-p_infin))/((initial_D**2)*RHO_P*RV*TEMPERATURE) #evaporation rate
    #beta = (8*molec_diff*(p_sat-p_infin))/((droplet_diameter(time-0.005, D_0)**2)*RHO_P*RV*TEMPERATURE) #using previous droplet size

    d_min = 0.44*initial_D
    d = max(initial_D*math.sqrt(max(1-beta*time,0)), d_min)

    return d

def terminal_velocity(time,initial_D):
    ''' This function estimates the terminal velocity in m/s of a respitory droplet as a function of
    the droplet's diamter "d" in micro meters. The terminal velocity also depends on defined constants and variables,
    gravity, drad coefficient, Reynolds number, and the density of the droplet and air.

    Parameters:
        d (float): A float representing the the diamter of the droplet which will be used to calculate the 
                Reynolds number and then used to find the drag coefficient "Cd"

    Returns: 
        v_t (float): v_t, a float value that represents the terminal velocity of the droplet, when Fdrag = Fgrav 
                    and neglects the Buoyant force. 

    '''
    if time <= 0:
        return (RHO_P*initial_D**2*G)/(18*math.pi*VISCOSITY) #Stoke's Law for small velocities

    #v_temp = terminal_velocity(time-0.005, initial_D)
    #d_before = droplet_diameter(time-0.005, initial_D)
    #reynolds_p = RHO_A*v_temp*d_before/VISCOSITY #reynolds number should be using diameter and v_t at previous iteration

    reynolds_p = RHO_A*V_X*d/VISCOSITY #reynolds number calculation
    drag_coef = 24*(1+0.15*reynolds_p**0.687)/reynolds_p #Drag Coefficient
    d = droplet_diameter(time,initial_D)
    v_t = math.sqrt((4*d*(RHO_D - RHO_A)*G)/(3*RHO_A*drag_coef))
    return v_t


def position(time,initial_D): 
    ''' This function estimates the horizontal and vertical distance the droplet has travelled at an inputted time for a 
    certain terminal velocity, horizontal velocity, initital X_0 and Z_0 

    Parameters:
        time (float): Time in seconds at which the x_d and z_d values are evaluated.

    Returns: 
        (x_d,z_d): a 2-tuple of float values containing the horizontal and vertical distance of the drop in meters
                   after t seconds.
    '''
    if time <= 0:
        return (X_0, Z_0)
        
    d = droplet_diameter(time,initial_D)
    v_t = terminal_velocity(d)

    time_to_hit_ground = Z_0/v_t
    x_d = X_0 + V_X*min(time, time_to_hit_ground)
    z_d = max(Z_0-v_t*time,0) #take into account droplet hitting the ground

    distance_tuple = (x_d,z_d)
    return distance_tuple

def concentration(time,x_away,initial_D): 
    ''' This function estimates the concentration of the droplets by modeling each breath as an ever expanding 
    Gaussian distribution.

    Parameters:
        time (float): Time in seconds at which the Gaussian distribution is evaluated.

    Returns: 
        conc_of_puff (float): a float value containing the concentration of the puff that interacts with a person X_AWAY 
        from the infected source.
    ''' 
    distance_tuple = position(time,initial_D)
    sigma = A*(x_away**B) #should this be x_d? not sure 
    x_d = distance_tuple[0]
    z_d = distance_tuple[1]
    conc_of_puff = (NUMBER_OF_DROPLETS/(math.sqrt(2*math.pi*sigma))**3)*math.exp((-1/(2*sigma**2))*((x_away-x_d)**2+z_d**2))
    return conc_of_puff

def exposure_per_breath(time,x_away,initial_D): 
    ''' This function estimates the dose of respiratory droplets that a person is exposed to by integrating the concentration 
    function above over the duration that the puff travels in the air. This function uses the python library scipy and the quad 
    function to evaluate the integral. 

    Parameters:
        time (float): Time in seconds that represents the upper limit at which the integral is evaluated.

    Returns: 
        exposure (2-tuple float): A 2-tuple of float value containing the concentration of the puff integrated over the inputted 
        time, and the error due to possible numerical error in the integrand from the use of quad 
    ''' 
    exposure = integrate.quad(concentration, 0, time, args=(x_away,initial_D,)) #keep x_away constant while integrating
    return exposure

def total_exposure(time,x_away=X_AWAY,initial_D=D_0):
    ''' This function estimates the total dosage a person is exposed to after many Guassian puffs by multiplying the
    respiratory rate x time to find the number of times the infected person has exhaled droplets into the air by exposure_per_breath
    function. 

    Parameters:
        time (float): Time in seconds for which the the total exposure is calculated.

    Returns: 
        total_dosage (float): a float value representing the total dosage a person is exposed to after several breaths are
        taken from an infected source. 
    ''' 
    exposure_tuple = exposure_per_breath(time,x_away,initial_D)
    number_of_breaths = RESPIRATORY_RATE*time
    total_dosage = exposure_tuple[0]*number_of_breaths
    #print(total_dosage)
    return total_dosage
    

if __name__ == '__main__':
    total_exposure(50,4)

    t = 3 #three iterations of time because only 3 initial_D's ?
    initial_D_list = list(np.arange(10*10**-6, 100*10**-6, 10**-6))

    x_away = [0.25,0.5,1,2,3]
    for x in x_away:
        exposure_array = []
        for init_D in initial_D_list:
            exposure = total_exposure(t,x,init_D)
            exposure_array.append(exposure)
        plt.plot(initial_D_list,exposure_array, label = "x_away = " + str(x))
    plt.xlabel('Droplet Size')
    plt.ylabel('Concentration of Droplets')
    plt.title('Concentration vs Droplet Size Graph')
    plt.legend()
    plt.show()
