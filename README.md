# aLIGO-wxPython
Graphical User Interface for Gravitational Wave Data Analysis. 

For an specific compact binary system of the same mass and at a given distance it calculates the Signal-to-Noise Ratio, Horizon Distance. Average Range, Time to innermost stable circular orbits (ISCO), and Time in view for three different detectors configuration (LIGO Hanford, LIGO Livingston and a Network of the two detectors). Three design parameters can be change: 1-Total Suspension Length, 2- Section Suspension Length and 3-Total mass of the suspension. Can select from 9 different masses, and 3 lengths for both the total and section length. Four source parameters can be changed: 1-Mass of the two binaries, 2-Distance of the binaries, 3-Position in the Sky (converted to earth longitude and latitude), 4-Maximun gravitaional wave (GW) frequency (GW freq. is twice the orbital frequency). It plots and calculate the amplitude spectral density (ASD) from the Seismic, Suspension Thermal, Coating Brownian and Quantum Noise. 


Assumptions
======= 
The code assumes polarization angle of 0.343, observation time at t0=9e8 (gps sec), and inclination of 0 (best scenario). It is a stationary phase approximation.

Prerequisite packages
=======
<!---
<dt>**LSC Algorithm Library Suite:**</dt>
* [LALSuite](https://www.lsc-group.phys.uwm.edu/daswg/projects/lalsuite.html). Full instruction to install can be found in (https://www.lsc-group.phys.uwm.edu/daswg/docs/howto/lal-install.html).
-->

<dt>**Python Packages:**</dt>
* [NumPy](http://www.numpy.org/)
* [Matplotlib](http://matplotlib.org/)
* [Basemap](http://matplotlib.org/basemap/)
* [wxPython](http://www.wxpython.org/)
* [Astropy](http://www.astropy.org/)


Usage:
=======
To execute clone the git repository and run the main program aLIGO:
* `python aLIGO.py`

Gallery
=======

**Screen Shot:**

![Alt text](https://cloud.githubusercontent.com/assets/8272801/9312808/16c1d406-4517-11e5-9607-ca9b75e49d61.png)

Future Work
=======
* Include gravity gradient noise
* Include time as a parameters
* More datapoints
* Include signal after ISCO

