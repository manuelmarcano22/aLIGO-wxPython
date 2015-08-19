#!/usr/bin/env python
import lal
import lalsimulation
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from astropy.time import Time


def antenna_response( t0, ra, dec, psi, det ):
    gps = lal.LIGOTimeGPS( t0 )
    gmst_rad = lal.GreenwichMeanSiderealTime(gps)

    # create detector-name map
    detMap = {'H1': lal.LALDetectorIndexLHODIFF, \
    'H2': lal.LALDetectorIndexLHODIFF, \
    'L1': lal.LALDetectorIndexLLODIFF, \
    'G1': lal.LALDetectorIndexGEO600DIFF, \
    'V1': lal.LALDetectorIndexVIRGODIFF, \
    'T1': lal.LALDetectorIndexTAMA300DIFF}

    try:
        detector=detMap[det]
    except KeyError:
        raise ValueError, "ERROR. Key %s is not a valid detector name." % (det)

    # get detector
    detval = lal.CachedDetectors[detector]

    response = detval.response

    # actual computation of antenna factors
    fp, fc = lal.ComputeDetAMResponse(response, ra, dec, psi, gmst_rad)
    ft = np.sqrt(fp**2 + fc**2)
    return ft
    
t0=900001970 #gps seconds
det = 'L1'

#Hanford
if det =='H1':
    Lambda = np.pi * 46.45 / 180; # latitude
    gamma = np.pi * 171.8 / 180; # orientation of arms
    lde = np.pi * 119.41 / 180;  # longitude 
    xi = np.pi / 2; # angle between arms
#Livinsgton
if det == 'L1':
    Lambda = np.pi * 30.56 / 180; # latitude
    gamma = np.pi * 243.0 / 180; # orientation of arms
    lde = np.pi * 90.77 / 180; # longitude
    xi = np.pi / 2; # angle between arms


psi = 0#0.343 #(radians) Polarization angle
t = Time(t0, format='gps')
gmst = t.sidereal_time('mean','greenwich').rad
lst = t.sidereal_time('mean',longitude=str(360-np.degrees(lde))+'d')
print 'in ', t.isot
print 'max at (ra, dec), delta:', 180-lst.deg, np.degrees(Lambda), ((180-lst.deg)-np.degrees(lde) )

Xlista = np.arange(0,360.,1.)
Ylista = np.arange(-90.,90.,1.)
Z = np.array([])


for i in Ylista:
    for j in Xlista:
        ztemp = antenna_response(t0, np.radians(j),np.radians(i), psi, det ) 
        Z = np.append(Z,ztemp)

X = np.arange(-180,180.,1.)
#X = np.arange(0,360.,1.)
Y = np.arange(-90.,90.,1.)
X, Y = np.meshgrid(X,Y)
#Z = np.reshape(Z, X.shape)
#
## create figure, axes instances.
fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])

## create Basemap instance.
m = Basemap(projection='moll',lon_0=0, lat_0=0)
#X, Y = m.makegrid(X,Y)
Z = np.reshape(Z, X.shape)

#X,Y= m(X,Y)
m.drawcoastlines()
im1 = m.pcolormesh(X,Y,Z,shading='flat',cmap=plt.cm.jet,latlon=True)
plt.show()
#
#
#
