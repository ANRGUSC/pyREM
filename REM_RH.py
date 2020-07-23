import math
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root
from scipy.integrate import simps

#np.seterr(all='warn')
np.seterr(all='raise')
gvt  = 0 

RHO_A = 1.21 #density of air in kg/m^3 
RHO_D = 1000 #density of droplet in kg/m^3
RHO = RHO_A
RHO_P = RHO_D
G = 9.81 #gravitational acceleration in m/s^2
VISCOSITY = 1.81*10**-5 #viscosity of air in Pa s
RV = 461.52 #J/kgK specific gas constant for water
D_0 = 9.3*10**-5 #initial diameter of droplet
A = 0.06 #given constant in dispersion coefficient equation
B = 0.92 #given constant in dispersion coefficient equation
NUMBER_OF_DROPLETS = 1 #number of droplets emitted (q)
X_0 = 0 #initial X position
Z_0 = 0 #Initial vertical position 
RESPIRATORY_RATE = 0.25 #breaths/second from avg of 15 bpm
RELATIVE_HUMMIDITY = 60 #relative hummidity
TEMPERATURE = 293.15 # ambient temperature in Kelvin
V_X = 1 #horizontal velocity 1m/s
X_AWAY = 2 #a distance X meters away from source 


def diameter_polynomial(time,r_h,initial_D):
    ''' This function estimates the droplet's diameter in micrometers as a function of time (s) that also depends
    on the relative hummidity (RH) and the temperature (T) in Kelvin which are set as constants in this program.
    Parameters:
        time (float): This parameter represents the time at which the diameter of the droplet will be calclulated, 
                    since the size of the droplet evaporates over time.
    Returns:
        d (float): Returns d, a float value representing the diameter of the droplet after it has 
                evaporated after t seconds. 
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
    global gvt
    if time <= 0:
        v_t = (RHO_P*initial_D**2*G)/(18*math.pi*VISCOSITY) #Stoke's Law for small velocities
    else:
        d = diameter_polynomial(time,r_h,initial_D) 
        n = 10.8*VISCOSITY*((RHO_A*d)/VISCOSITY)**0.687 
        p = 4*(d**2)*(RHO_D-RHO_A)*G 
        m = 72*VISCOSITY
        try: 
            roots = root(lambda v: n*v**(2.687)+m*v**2-p*v,0.1)
            v_t = roots.x[0]
            gvt = v_t
        except: 
            v_t = gvt 
    #print(v_t)
    return v_t


def position(time,r_h,initial_D): 
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
    ''' This function estimates the concentration of the droplets by modeling each breath as an ever expanding 
    Gaussian distribution.
    Parameters:
        time (float): Time in seconds at which the Gaussian distribution is evaluated.
    Returns: 
        conc_of_puff (float): a float value containing the concentration of the puff that interacts with a person X_AWAY 
        from the infected source.
    ''' 
    distance_tuple = position(time,r_h,initial_D)
    x_d = distance_tuple[0]
    z_d = distance_tuple[1]
    sigma = A*(x_d**B) 
    conc_of_puff = (NUMBER_OF_DROPLETS/((math.sqrt(2*math.pi)*sigma))**3)*math.exp((-1/(2*sigma**2))*((X_AWAY-x_d)**2+z_d**2))

    #print(time, conc_of_puff)

    return conc_of_puff

def exposure_per_breath(time,r_h,initial_D): 
    ''' This function estimates the dose of respiratory droplets that a person is exposed to by integrating the concentration 
    function above over the duration that the puff travels in the air. This function uses the python library scipy and the quad 
    function to evaluate the integral. 
    Parameters:
        time (float): Time in seconds that represents the upper limit at which the integral is evaluated.
    Returns: 
        exposure (2-tuple float): A 2-tuple of float value containing the concentration of the puff integrated over the inputted 
        time, and the error due to possible numerical error in the integrand from the use of quad 
    ''' 
    #exposure = integrate.quad(concentration, 0, time, args=(r_h,initial_D,), limit=100) #keep x_away constant while integrating
    #return exposure

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
    #print(exposure)

    return exposure 

    
def total_exposure(time,r_h=RELATIVE_HUMMIDITY,initial_D=D_0):
    ''' This function estimates the total dosage a person is exposed to after many Guassian puffs by multiplying the
    respiratory rate x time to find the number of times the infected person has exhaled droplets into the air by exposure_per_breath
    function. 
    Parameters:
        time (float): Time in seconds for which the the total exposure is calculated.
    Returns: 
        total_dosage (float): a float value representing the total dosage a person is exposed to after several breaths are
        taken from an infected source. 
    ''' 
    exposure_tuple = exposure_per_breath(time,r_h,initial_D)
    number_of_breaths = RESPIRATORY_RATE*time
    #total_dosage = exposure_tuple[0]*number_of_breaths
    total_dosage = exposure_tuple[0]
    #print(total_dosage)

    return total_dosage
    

if __name__ == '__main__':
    gvt = 0
    #total_exposure(5)
    #terminal_velocity(5,60,D_0)

    time_array = []
    z_d_array = []
    vt_array = []

    #for t in range(0,5):
    for i in range(1,50):
        t = (i)/(10.0)
        d = 93*10**-6
        distance_tuple = position(t,60,d)
        vt = terminal_velocity(t,60,d)
        z_d = distance_tuple[1]
        z_d_array.append(z_d)
        time_array.append(t)
        vt_array.append(vt)
       
    #plt.plot(time_array,z_d_array)
    #plt.xlabel('Time')
    #plt.ylabel('Z position')
    #plt.title('Z pos vs Time @93um')
    #plt.show()

    plt.plot(time_array,vt_array)
    plt.xlabel('Time')
    plt.ylabel('Terminal Velocity')
    plt.title('Velocity vs Time @93um')
    plt.show()

'''
# plot for trajectories
    initial_D_array = [85,89,93,95,97]

    for init_D in initial_D_array:
        x_d_array = []
        z_d_array = []
        for t in range(0,10):
            d = init_D*10**-6
            distance_tuple = position(t,60,d)
            x_d = distance_tuple[0]
            z_d = distance_tuple[1]
            x_d_array.append(x_d)
            z_d_array.append(z_d)
        #print(x_d_array)
        #print(z_d_array)
        plt.plot(x_d_array,z_d_array, label = "D_0 = " + str(init_D))
    plt.xlabel('X position')
    plt.ylabel('Z position')
    plt.title('Trajectories')
    plt.legend()
    plt.show()
'''

''' #plot for conc vs humidity
    t = 5
    initial_D_list = list(np.arange(1*10**-6, 1*10**-4, 2*10**-6))
    #hummidity = [50,60,70,80]
    hummidity = [60]
    for r in hummidity:
        exposure_array = []
        for init_D in initial_D_list:
            print(init_D)
            exposure = exposure_per_breath(t,r,init_D)
            print(exposure)
            exposure_array.append(exposure)
        plt.plot(initial_D_list,exposure_array, label = "RH = " + str(r))
    plt.xlabel('Droplet Size')
    plt.ylabel('Concentration of Droplets')
    plt.title('Concentration vs Droplet Size Graph')
    plt.legend()
    plt.show() 
'''

