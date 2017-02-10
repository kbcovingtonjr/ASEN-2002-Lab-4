Wind Tunnel Lab
------

After taking data on an airfoil in a low-speed subsonic wind tunnel here at CU Boulder, I wrote this code to analyze my data in conjunction with all of the data taken by my class.

The goal was to find the coefficients of lift and drag for a Clark Y-14 airfoil.
![ClCdGraphs](/Figures/ClCdGraphs.jpg?raw=true "ClCdGraphs")


This is compared to the official NACA data for the airfoil.
![NACAcomparison](/Figures/NACA_comparison.jpg?raw=true "NACA data comparison")


In addition, I plotted the pressure distribution across the airfoil.
![PressureDistribution](/Figures/pressure_distribution.jpg?raw=true "Pressure distribution")

More In-Depth
------
Data was taken by many groups of students within the CU Boulder Aerospace Engineering Class of 2019. These data, together, constitute comprehensive data such on the low-speed aerodynamics of the NACA Clark Y-14 airfoil. Data was taken at wind speeds of 10, 20, and 30 m/s and angles of attack ranging from -15 to 15 degrees.

This code:
* Combines, averages, and parses wind tunnel data from several groups' .csv files.
* Calculates coefficients of pressure for each pressure transducer on the airfoil.
* Uses these Cp values to numerically approximate the pressure distribution over the airfoil and calculates coefficients of lift and drag.
* Plots pressure distribution curves.
* Reports findings in an output .txt file.
