import matplotlib.pyplot as plt
import numpy as np 
from REM_temp import diameter_polynomial,terminal_velocity,position,concentration,exposure_per_breath,total_exposure


if __name__ == '__main__':

#   plot for concentration vs droplet size with temp varied
    t = 5
    initial_D_list = list(np.arange(1*10**-6, 1*10**-4, 2*10**-6))
    temperature = [0,10,20,30]
    for temp in temperature:
        exposure_array = []
        temp_k = temp+273.15
#        print(temp_k)
        for init_D in initial_D_list:
#            print(init_D)
            exposure = total_exposure(t,temp_k,init_D)/(12.3487382) #normalize the curves
            exposure_array.append(exposure)
        plt.plot(initial_D_list,exposure_array, label = "T = " + str(temp))
    plt.xlabel('Droplet Size')
    plt.ylabel('Concentration of Droplets')
    plt.title('Concentration vs Droplet Size Graph')
    plt.legend()
    plt.show() 
