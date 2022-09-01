# Coccolith_CIVEs
Codes for exploring Coccolith Carbon Isotope Vital Effects (CIVEs)

This respository contains Matlab codes developed for use in the following papers to model carbon isotope vital effects (CIVEs) in coccolith calcite :

McClelland et al., 2017, Nature Communications.   
link: https://www.nature.com/articles/ncomms14511  
Questions to:  hmcclelland@unimelb.edu.au

Claxton et al., 2022, Nature Geoscience.  
link: https://www.nature.com/articles/s41561-022-01006-0  
Questions to: louis.m.claxton@gmail.com OR hmcclelland@unimelb.edu.au

# Set-up
Master_file.m is the file to run in Matlab. You may need to make sure Matlab is up to date (2021), and install the following toolboxes: Statistics and Machine Learning Toolbox; and Parallel Computing Toolbox. 

# Files

- Master_File.m - runs optimisation routine presented in Claxton et al. 2022.

- functions/
	- cocco_iso_concs.m - Analytical solution to system of ODEs (McClelland at al. 2017)
	- ISOMOD_analytical_fast.m - Matlab implementation of isotope flux model (McClelland at al. 2017). Vectorized for speed.
	- dDICs_fast.m - Fast calculation of equilibrium carbonate chemistry parameters (concs and isotopes)
	- forward_13C.m - Calculates the isotopic composition of a multi-species size fraction, minus the mean across all size fractions.  
	- forward_13C_absolute.m - As forward_13C.m but does not minus the mean for comparison with the foram record.  
	- minfunCO2PP - minimization function called by Master_File.m to iteratively constrain CO2 and PP. 
