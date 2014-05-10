spuc_filter
===========

A Qt5 based GUI demonstration of the various filtering capabilities of the SPUC library

This uses 

QCustomPlot library (www.qcustomplot.com) - GPL v3

Requirements:

		boost
		cmake
		qt5
		spuc

		(see .travis.yml for install steps if needed)
		cd build
		cmake ..
		make


Use the mouse to select cut-off, etc

Simulataneously enable multiple filters for comparison

Filter types in SPUC (not all in this app)
	* Butterworth
	* Chebyshev
	* Elliptic
	* Maximally flat FIR
	* Remez Equiripple
	* Raised Cosine FIR
	* Gaussian FIR
	* CIC
	* Notch filter
	* Cut/Boost Filter
	* Halfband/Subband IIR filters
	   
![Demo App](app.png "Demo App")
