import matplotlib.pyplot as plt
import numpy as np 
from REM_Xaway import diameter_polynomial,terminal_velocity,position,concentration,exposure_per_breath,total_exposure
 

if __name__ == '__main__':

#   plot for concentration vs droplet size with x_away varied
    t = 5
    initial_D_list = list(np.arange(1*10**-6, 1*10**-4, 2*10**-6))
    x_away = [0.25,0.5,1,2]
    for x in x_away:
        exposure_array = []
#        print(x)
        for init_D in initial_D_list:
#            print(init_D)
            exposure = total_exposure(t,x,init_D)/(15.43592275) #normalize the curves
            exposure_array.append(exposure)
        plt.plot(initial_D_list,exposure_array, label = "x_away = " + str(x))
    plt.xlabel('Droplet Size')
    plt.ylabel('Concentration of Droplets')
    plt.title('Concentration vs Droplet Size Graph')
    plt.legend()
    plt.show() 
