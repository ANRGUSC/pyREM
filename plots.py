import matplotlib.pyplot as plt
import numpy as np 
from REM import diameter_polynomial,terminal_velocity,position,concentration,exposure_per_breath,total_exposure

RELATIVE_HUMIDITY = 60 #default relative humidity
TEMPERATURE = 293.15 #default ambient temperature in Kelvin
X_AWAY = 2 #default distance 2 meters away from source 
N_FACTOR = total_exposure(5,X_AWAY,TEMPERATURE,RELATIVE_HUMIDITY,1*10**-6) #dose of smaller droplets used to normalize curves

def proximity_plot(time):
    initial_D_list = list(np.arange(1*10**-6, 1*10**-4, 2*10**-6))
    x_away = [0.25,0.5,1,2]
    for x in x_away:
        exposure_array = []
        print(x)
        for init_D in initial_D_list:
            print(init_D)
            exposure = total_exposure(time,x,TEMPERATURE,RELATIVE_HUMIDITY,init_D)/N_FACTOR #normalize the curves
            exposure_array.append(exposure)
        plt.plot(initial_D_list,exposure_array, label = "x_away = " + str(x))
    plt.xlabel('Droplet Size')
    plt.ylabel('Concentration of Droplets')
    plt.title('Concentration vs Droplet Size Graph')
    plt.legend()
    plt.show() 

    return

def temp_plot(time):
    initial_D_list = list(np.arange(1*10**-6, 1*10**-4, 2*10**-6))
    temperature = [0,10,20,30]
    for temp in temperature:
        exposure_array = []
        temp_k = temp+273.15
        print(temp_k)
        for init_D in initial_D_list:
            print(init_D)
            exposure = total_exposure(time,X_AWAY,temp_k,RELATIVE_HUMIDITY,init_D)/N_FACTOR #normalize the curves
            exposure_array.append(exposure)
        plt.plot(initial_D_list,exposure_array, label = "T = " + str(temp))
    plt.xlabel('Droplet Size')
    plt.ylabel('Concentration of Droplets')
    plt.title('Concentration vs Droplet Size Graph')
    plt.legend()
    plt.show() 

    return 

def humidity_plot(time):
    initial_D_list = list(np.arange(1*10**-6, 1*10**-4, 2*10**-6))
    hummidity = [50,60,70,80]
    for r in hummidity:
        exposure_array = []
        print(r)
        for init_D in initial_D_list:
            print(init_D)
            exposure = total_exposure(time,X_AWAY,TEMPERATURE,r,init_D)/N_FACTOR #normalize the curves
            exposure_array.append(exposure)
        plt.plot(initial_D_list,exposure_array, label = "RH = " + str(r))
    plt.xlabel('Droplet Size')
    plt.ylabel('Concentration of Droplets')
    plt.title('Concentration vs Droplet Size Graph')
    plt.legend()
    plt.show() 

    return 

#functions for intermediary plots
def diameter_plot(time=50):
    initial_D_array = [10,33,75,92]

    for init_D in initial_D_array:
        time_array = []
        d_array = []
        for i in range(1,time):
            t = (i)/(10.0)
            time_array.append(t)
            d = init_D*10**-6
            dm = diameter_polynomial(t,TEMPERATURE,RELATIVE_HUMIDITY,d)
            d_array.append(dm)
        plt.plot(time_array,d_array, label = "D_0 = " + str(init_D))
    plt.xlabel('Time')
    plt.ylabel('Diameter (m)')
    plt.title('Diameter vs Time')
    plt.legend()
    plt.show()

    return 

def velocity_plot(time):
    initial_D_array = [10,33,75,92]

    for init_D in initial_D_array:
        time_array = []
        vt_array = []
        for i in range(1,time):
            t = (i)/(10.0)
            time_array.append(t)
            d = init_D*10**-6
            vt = terminal_velocity(t,TEMPERATURE,RELATIVE_HUMIDITY,d)
            vt_array.append(vt)
        plt.plot(time_array,vt_array, label = "D_0 = " + str(init_D))                                   
    plt.xlabel('Time')
    plt.ylabel('Terminal Velocity')
    plt.title('Velocity vs Time')
    plt.legend()
    plt.show()

    return

def x_pos_plot(time):
    time_array = []
    x_d_array =[]
    for i in range(1,time):
        t = (i)/(10.0)
        time_array.append(t)
        d = 85*10**-6
        distance_tuple = position(t,TEMPERATURE,RELATIVE_HUMIDITY,d)
        x_d = distance_tuple[0]
        x_d_array.append(x_d)
    plt.plot(time_array,x_d_array, label = "x_d vs time for all D_0's")
    plt.xlabel('Time')
    plt.ylabel('X position')
    plt.title('X Position vs Time')
    plt.legend()
    plt.show()

    return

def z_pos_plot(time):
    initial_D_array = [10,33,75,92]

    for init_D in initial_D_array:
        time_array = []
        z_d_array = []
        for i in range(1,time):
            t = (i)/(10.0)
            time_array.append(t)
            d = init_D*10**-6
            distance_tuple = position(t,TEMPERATURE,RELATIVE_HUMIDITY,d)
            z_d = distance_tuple[1]
            z_d_array.append(z_d)
        plt.plot(time_array,z_d_array, label = "D_0 = " + str(init_D)) 
    plt.xlabel('Time')
    plt.ylabel('Z position')
    plt.title('Z Position vs Time')
    plt.legend()
    plt.show()
    return

def trajectories(time=10):
    initial_D_array = [10,33,75,92]

    for init_D in initial_D_array:
        x_d_array = []
        z_d_array = []
        for t in range(0,time):
            d = init_D*10**-6
            distance_tuple = position(t,TEMPERATURE,RELATIVE_HUMIDITY,d)
            x_d = distance_tuple[0]
            z_d = distance_tuple[1]
            x_d_array.append(x_d)
            z_d_array.append(z_d)
        plt.plot(x_d_array,z_d_array, label = "D_0 = " + str(init_D))
    plt.xlabel('X position')
    plt.ylabel('Z position')
    plt.title('Trajectories')
    plt.legend()
    plt.show()

    return

def concentration_plot(time):
    initial_D_array = [10,33,75,92]

    for init_D in initial_D_array:
        time_array = []
        conc_array = []
        for i in range(1,time):
            t = (i)/(10.0)
            time_array.append(t)
            d = init_D*10**-6
            conc = concentration(t,X_AWAY,TEMPERATURE,RELATIVE_HUMIDITY,d)
            conc_array.append(conc)
        plt.plot(time_array,conc_array, label = "D_0 = " + str(init_D))
    plt.xlabel('Time')
    plt.ylabel('Concentration')
    plt.title('Concentration vs Time')
    plt.legend()
    plt.show()

    return

#example usage, for testing
if __name__ == '__main__':
#    proximity_plot(5)
#    temp_plot(5)
#    humidity_plot(5)

#intermediary plots: D(t), v(t), z(t), x(t), and trajectory
#    diameter_plot(50)
#     velocity_plot(50)
#     x_pos_plot(50)
#     z_pos_plot(50)
#    trajectories(10)
#    concentration_plot(50)
