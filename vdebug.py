import matplotlib.pyplot as plt
import numpy as np 
from scipy.optimize import root
import math,pdb

RHO_A = 1.21 #density of air in kg/m^3 
RHO_D = 1000 #density of droplet in kg/m^3
RHO = RHO_A
RHO_P = RHO_D
G = 9.81 #gravitational acceleration in m/s^2
VISCOSITY = 1.81*10**-5 #viscosity of air in Pa s
RV = 461.52 #J/kgK specific gas constant for water
D_0 = 2.45*10**-5 #initial diameter of droplet in micro meters
RELATIVE_HUMIDITY = 60 #default relative humidity
TEMPERATURE = 293.15 #default ambient temperature in Kelvin

molec_diff = (2.16*10**-5)*(TEMPERATURE/273.15)**1.8 # molecular diffusivity of water vapor
p_sat = 611.21*math.exp((19.843-(TEMPERATURE/234.5))*((TEMPERATURE-273.15)/(TEMPERATURE-16.01))) # saturation water vapor pressure
p_infin = p_sat*RELATIVE_HUMIDITY/100 # ambient water vapor pressure
t_crit = (RHO_P*RV*TEMPERATURE*(D_0**2))/(32*molec_diff*(p_sat-p_infin))

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
    print(d,t_crit)
    return d

def terminal_velocity(time,temp,r_h,initial_D):
    if time <= 0:
        v_t = (RHO_P*initial_D**2*G)/(18*math.pi*VISCOSITY) #Stoke's Law for small velocities
    else:
        d = diameter_polynomial(time,temp,r_h,initial_D) 
        n = 10.8*VISCOSITY*((RHO_A*d)/VISCOSITY)**0.687 
        p = 4*(d**2)*(RHO_D-RHO_A)*G 
        m = 72*VISCOSITY
        roots = root(lambda v: n*v**(2.687)+m*v**2-p*v,1)
        v_t = roots.x[0]
        print(v_t)

    return v_t

if __name__ == '__main__':
    diameter_polynomial(t_crit,293.15,60,D_0)

'''  
    #initial_D_array = [78,85,90,96]

    for init_D in initial_D_array:
        time_array = []
        vt_array = []
        for i in range(1,50):
            t = (i)/(10.0)
            time_array.append(t)
            d = init_D*10**-6
            print(t)
            vt = terminal_velocity(t,TEMPERATURE,RELATIVE_HUMIDITY,d) #BUGGGG
            vt_array.append(vt)
        plt.plot(time_array,vt_array, label = "D_0 = " + str(init_D))                                   
    plt.xlabel('Time')
    plt.ylabel('Terminal Velocity')
    plt.title('Velocity vs Time')
    plt.legend()
    plt.show()

   '''
