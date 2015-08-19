# aLIGO-wxPython
**Graphical User Interface for Gravitational Wave Data Analysis.**

For an specific compact binary system of the same mass and at a given distance it calculates the Signal-to-Noise Ratio, Horizon Distance. Average Range, Time to innermost stable circular orbits (ISCO), and Time in view for three different detectors configuration (LIGO Hanford, LIGO Livingston and a Network of the two detectors). For the fist prototype of the GUI we will only consider the effect of changing the suspension design parameters. Advanced LIGO has a quadruple suspension and we can change three important things on the Suspension:

* Total mass of the suspension. 9 possible masses from 40 kg to 120 kg
* Total length. 3 possible configurations (1.6, 1.87 and 2.14 meters).
* Sections length (ribbons). Also 3 possibilities (0.6, 0.85 and 1.1 meters).

Four source parameters can also be changed: 

* Mass of the two binaries
* Distance of the binaries
* Position of the source (converted to Earth-fixed coordinates)
* Maximun Gravitational Wave frequency

It plots and calculate the effective induced strain and amplitude spectral density (ASD) from the Seismic, Suspension Thermal, Coating Brownian and Quantum Noise.


Assumptions
======= 
The code assumes polarization angle of 0., observation time at t0=9e8 (gps sec), and inclination of 0 (best scenario). It is a stationary phase approximation and assumes no spinning binaries in a circular orbit. 

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

Usage:
=======
To execute clone the git repository and run the main program aLIGO:
* `python aLIGO.py`


Report:
=======
Find a short report on the project on the report folder. **Warning: It is a work in progress and the report is not complete**


Gallery
=======

**Screen Shot:**

![Alt text](https://cloud.githubusercontent.com/assets/8272801/9312808/16c1d406-4517-11e5-9607-ca9b75e49d61.png)

Future Work
=======
* Include gravity gradient noise
* Include time as a parameters
* More datapoints for the sliders
* Include contribution to SNR beyond ISCO
* Inclination as parameter or randomize
* Include eccentric and spinning binaries

