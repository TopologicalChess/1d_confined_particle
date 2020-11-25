This program assumes that gnuplot 5.2 is up to be called by the system!
Author: Victor Gustavo May Custodio
Date:  24/11/2020

1D confined particle----------------------------------------------------
	This program is intended to solve the Schrödinger equation for a given
	potential V(x) defined on some interval I=[a,b]. We assume the particle
	is confined in this region of space by the boundary conditions 
	Psi(a)=Psi(b)=0.
	
	There are some example potential included, but one could add a .dat file
	to this folder and the program would ask if what potential to be used
	The introduced potential should be formatted in the following way:

		1.- two columns (x,v(x))
		2.- the first column should be in meters
		3.- the second column should be the measured potential
		    in Joules divided by the factor hbar^2/2m
		
	The program then outputs the wavenumbers.dat, eigenstates.dat and plotting.dat
	which contain the information of the simulation. If gnuplot is running, then the
	program calls gnuplot for plotting the first 4 states of the 1-D system. The settings
	of this final plotting can be adjusted in the plots.plt file. The final file is the 
	firststates.png which contains the plotted image from gnuplot.

	Gnuplot script "plots.plt" is adjusted to run Quantum Harmonic Oscillator from -10.0 to 10.0.
	 