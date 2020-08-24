# pyREM

This repository contains open source python code that can be used to quantify an individual’s accumulated exposure to respiratory droplets containing SARS-CoV-2 viral RNA over a given period of time. Our code is based on Professor Gavin Buxton’s “Spreadsheet Model of COVID-19 Transmission” [1] and is a python based implementation of this paper. It contains the following programs.
 
- REM.py: This python file contains a function `total_exposure()` which takes time, distance, temperature and humidity as inputs and outputs the accumulated viral exposure after t seconds. If certain input parameters are not given, the function uses default values specified in the original paper: 
`x_away=2,temp=273.15,r_h=60, initial_D=1*10**-5`

Example usage: `total_exposure(5)`

- plots.py: This python file imports the functions from REM.py and contains a series of functions that will generate the figures of our numerical results. Each function takes in time as an input and outputs the corresponding plot for the respective variable(s). 

Example usages:  
`proximity_plot(5),`
`temp_plot(5),`
`humidity_plot(5),`   
`diameter_plot(50),`
`velocity_plot(50),`
`x_pos_plot(50),`
`z_pos_plot(50),`
`trajectories(10),`
`concentration_plot(50)`
 
 A more detailed description of PyREM can be found in [2].
 
 
Contributors: Radhika Bhuckory (bhuckory@usc.edu) and Bhaskar Krishnamachari (bkrishna@usc.edu), University of Southern California


[1] Buxton, Gavin. Spreadsheet Model of COVID-19 Transmission: Evaporation and Dispersion of Respiratory Droplets, Spreadsheets in Education, vol 12, no. 2, May 2020. Available online at: https://sie.scholasticahq.com/article/12861-spreadsheet-model-of-covid-19-transmission-evaporation-and-dispersion-of-respiratory-droplets

[2] Radhika Bhuckory, Bhaskar Krishnamachari, PyREM: Implementing a Computational Model of Airborne Respiratory Droplet-based Virus Transmission, August 2020. Available at: https://github.com/ANRGUSC/pyREM/blob/master/PyREM_TechReport_August2020.pdf



