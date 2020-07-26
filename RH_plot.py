import matplotlib.pyplot as plt
import numpy as np 
from REM_RH import diameter_polynomial,terminal_velocity,position,concentration,exposure_per_breath,total_exposure


if __name__ == '__main__':

#   plot for concentration vs droplet size with humidity varied
    t = 5
    initial_D_list = list(np.arange(1*10**-6, 1*10**-4, 1*10**-5))
    hummidity = [50,60,70,80]
    for r in hummidity:
        exposure_array = []
        print(r)
        for init_D in initial_D_list:
            print(init_D)
            exposure = total_exposure(t,r,init_D)/(15.435922858566114)
            exposure_array.append(exposure)
        plt.plot(initial_D_list,exposure_array, label = "RH = " + str(r))
    plt.xlabel('Droplet Size')
    plt.ylabel('Concentration of Droplets')
    plt.title('Concentration vs Droplet Size Graph')
    plt.legend()
    plt.show() 



