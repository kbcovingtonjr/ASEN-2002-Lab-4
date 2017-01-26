I wrote this code to analyze wind tunnel data taken by multiple groups.
These data, together, constitute comprehensive data such on the
low-speed aerodynamics of the NACA Clark Y-14 airfoil. Data was taken
at wind speeds of 10, 20, and 30 m/s and angles of attack ranging from
-15 to 15 degrees.

This code performs the following:
	-Combines, averages, and parses wind tunnel data from several 
	 groups' .csv files.
	-Calculates coefficients of pressure for each pressure 
	 transducer on the airfoil.
	-Uses these Cp values to numerically approximate the pressure
	 distribution over the airfoil and calculates coefficients of 
	 lift and drag.
	-Plots pressure distribution curves.
	-Reports findings in an output .txt file.
