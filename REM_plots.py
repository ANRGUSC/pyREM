import matplotlib.pyplot as plt
import numpy as np 
from REM import diameter_polynomial,terminal_velocity,position,concentration,exposure_per_breath,total_exposure

if __name__ == '__main__':

# D(t), v(t), z(t), x(t), and trajectory plots
    time_array = []
    x_d_array =[]
    z_d_array = []
    vt_array = []
    d_array = []
  
    for i in range(1,50):
        t = (i)/(10.0)
        d = 93*10**-6
        time_array.append(t)

        distance_tuple = position(t,d)
        x_d = distance_tuple[0]
        x_d_array.append(x_d)

        z_d = distance_tuple[1]
        z_d_array.append(z_d)

        vt = terminal_velocity(t,d)
        vt_array.append(vt)

        dm = diameter_polynomial(t,d)
        d_array.append(dm)

# X position vs Time Plot
    #plt.plot(time_array,x_d_array, label = "D_0 = 93um")
    #plt.xlabel('Time')
    #plt.ylabel('X position')
    #plt.title('X Position vs Time')
    #plt.legend()
    #plt.show()

# Z position vs Time plot       
    #plt.plot(time_array,z_d_array, label = "D_0 = 93um,  RH = 60,  T = 20") 
    #plt.xlabel('Time')
    #plt.ylabel('Z position')
    #plt.title('Z Position vs Time')
    #plt.legend()
    #plt.show()

# Terminal Velocity vs Time plot
    #plt.plot(time_array,vt_array, label = "D_0 = 93um,  RH = 60,  T = 20")                                     
    #plt.xlabel('Time')
    #plt.ylabel('Terminal Velocity')
    #plt.title('Velocity vs Time')
    #plt.legend()
    #plt.show()

# Diameter vs Time plot
    #plt.plot(time_array,d_array, label = "D_0 = 93um,  RH = 60,  T = 20") 
    #plt.xlabel('Time')
    #plt.ylabel('Diameter (m)')
    #plt.title('Diameter vs Time')
    #plt.legend()
    #plt.show()

# plot for trajectories
    initial_D_array = [85,89,93,95,97]

    for init_D in initial_D_array:
        x_d_array = []
        z_d_array = []
        for t in range(0,10):
            d = init_D*10**-6
            distance_tuple = position(t,d)
            x_d = distance_tuple[0]
            z_d = distance_tuple[1]
            x_d_array.append(x_d)
            z_d_array.append(z_d)
#        plt.plot(x_d_array,z_d_array, label = "D_0 = " + str(init_D))
#    plt.xlabel('Z position')
#    plt.ylabel('Z position')
#    plt.title('Trajectories')
#    plt.legend()
#    plt.show()
