import matplotlib.pyplot as plt
import numpy as np 
from scipy.optimize import root
from REM import diameter_polynomial
import pdb

RHO_A = 1.21 #density of air in kg/m^3 
RHO_D = 1000 #density of droplet in kg/m^3
RHO = RHO_A
RHO_P = RHO_D
G = 9.81 #gravitational acceleration in m/s^2
VISCOSITY = 1.81*10**-5 #viscosity of air in Pa s
RV = 461.52 #J/kgK specific gas constant for water
D_0 = 1.0*10**-5 #initial diameter of droplet in micro meters
RELATIVE_HUMIDITY = 60 #default relative humidity
TEMPERATURE = 293.15 #default ambient temperature in Kelvin

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

    return v_t

if __name__ == '__main__':
    
    #initial_D_array = [78,85,90,96]
    initial_D_array = [78]

    for init_D in initial_D_array:
        time_array = []
        vt_array = []
        for i in range(1,50):
            t = (i)/(10.0)
            time_array.append(t)
            d = init_D*10**-6
            vt = terminal_velocity(t,TEMPERATURE,RELATIVE_HUMIDITY,d)
            print(vt)
#            vt_array.append(vt)
#        plt.plot(time_array,vt_array, label = "D_0 = " + str(init_D))                                   
#    plt.xlabel('Time')
#    plt.ylabel('Terminal Velocity')
#    plt.title('Velocity vs Time')
#    plt.legend()
#    plt.show()
