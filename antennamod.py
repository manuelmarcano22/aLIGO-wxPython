#!/usr/bin/env python

#from astropy.time import Time
import numpy as np

def antenna_pattern(det, ra, dec, psi, t):

    """Compute the antenna pattern for a given detector (site), 
    given a sky position in right ascension and declination (in rad), 
    and a polarisation angle psi. It outputs the pls and cross polarisation responses added in quadrature
    This code is adapted from the  matlab function by Dr. Matt Pitkin.
    """
    #Hanford
    if det ==0:
        Lambda = np.pi * 46.45 / 180; # latitude
        gamma = np.pi * 171.8 / 180; # orientation of arms
        lde = np.pi * 119.41 / 180;  # longitude 
        xi = np.pi / 2; # angle between arms
    #Livinsgton
    if det == 1:
        Lambda = np.pi * 30.56 / 180; # latitude
        gamma = np.pi * 243.0 / 180; # orientation of arms
        lde = np.pi * 90.77 / 180; # longitude
        xi = np.pi / 2; # angle between arms
        
    alpha = ra
    delta = dec

    #Get local sideral time
    #Convert gps time to LST
#    t = Time(t, format='gps')
#    gmst = t.sidereal_time('mean','greenwich').rad
#    lst = t.sidereal_time('mean',longitude=str(360-np.degrees(lde))+'d').rad
    lst = 1.5586
    #Convert to GMST from greewich:
    #lst = divmod(gmst-lde,2*np.pi)[1]
    gamma2 = 2 * gamma
    lambda2 = 2 * Lambda
    delta2 = 2 * delta
    s2g = np.sin(gamma2)
    c2g = np.cos(gamma2)
    s2l = np.sin(lambda2)
    sl = np.sin(Lambda)
    c2l = np.cos(lambda2)
    cl = np.cos(Lambda)
    s2d = np.sin(delta2)
    sd = np.sin(delta)
    c2d = np.cos(delta2)
    cd = np.cos(delta)
    alst = alpha - lst
    a = (1/16.) * s2g * (3 - c2l) * (3 - c2d) * \
            np.cos(2*(alst)) -\
            (1/4.) * c2g * sl * (3 - c2d) * \
            np.sin(2*(alst)) + \
            (1/4.) * s2g * s2l * s2d * \
            np.cos(alst) - \
            (1/2.) * c2g * cl * s2d * \
            np.sin(alst) + \
            (3/4.) * s2g * cl**2 * cd**2


    b = c2g * sl * sd * \
            np.cos(2*(alst)) + \
            (1/4.) * s2g * (3 - c2l) * sd * \
            np.sin(2*(alst)) + \
            c2g * cl * cd * \
            np.cos(alst) + \
            (1/2.) * s2g * s2l * cd * \
            np.sin(alst)
    
    Fp = np.sin(xi) * (a * np.cos(2*psi) + b * np.sin(2*psi))
    Fc = np.sin(xi) * (b* np.cos(2*psi) - a * np.sin(2*psi))
    Ftotal = np.hypot(Fp, Fc)

    return Ftotal
