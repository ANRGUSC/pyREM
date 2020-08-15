# pyREM

This repository contains open source python code that can be used to quantify an individual’s accumulated exposure to respiratory droplets containing SARS-CoV-2 viral RNA over a given period of time. Our code is based on Professor Gavin Buxton’s “Spreadsheet Model of COVID-19 Transmission” [1] and is a python based implementation of this paper. It contains the following programs.
 
- REM.py: This python file contains a function `total_exposure()` which takes time, distance, temperature and humidity as inputs and outputs the accumulated viral exposure after t seconds. If specific input parameters are not given, the function uses the default values defined in the original paper: 
`x_away=2,temp=273.15,r_h=60, initial_D= 10`

Example usage: `total_exposure(5)`

- plots.py: This python file imports the functions from REM.py and contains a series of functions that will generate figures of the numerical results when called. Each function takes in time as an input and outputs the corresponding plot for the respective variable(s). 

Example usages:  
`proximity_plot(5)
temp_plot(5)
humidity_plot(5)   
diameter_plot(50)
velocity_plot(50)
x_pos_plot(50)
z_pos_plot(50)
trajectories(10)
concentration_plot(50)`
 
 
Contributors: Radhika Bhuckory (bhuckory@usc.edu) and Bhaskar Krishnamachari (bkrishna@usc.edu), University of Southern California


[1] Buxton, Gavin, Spreadsheet Model of COVID-19 Transmission: Evaporation and Dispersion of Respiratory Droplets (April 22, 2020). Available at SSRN: https://ssrn.com/abstract=3582665 or http://dx.doi.org/10.2139/ssrn.3582665

