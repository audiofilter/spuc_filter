spuc_filter
===========

A Qt5 based GUI demonstration of the various filtering capabilities of the SPUC library

This uses 

QCustomPlot library (www.qcustomplot.com) - GPL v3

Requirements:

		cmake
		qt5
		spuc

		(see .travis.yml or .travis.linux for install steps if needed)
		cd build
		cmake ..
		make


Use the mouse to select cut-off, etc

Simulataneously enable multiple filters for comparison

Filter types in SPUC (not all in this app):

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
	   

### Build status - Mac Os X, clang - Automated Travis Build
[![Build Status](https://travis-ci.org/audiofilter/spuc_filter.png)](https://travis-ci.org/audiofilter/spuc_filter)

[![Gitter chat](https://badges.gitter.im/audiofilter/spuc.png)](https://gitter.im/audiofilter/spuc)

![Demo App](app.png "Demo App")
