#!/usr/bin/env python
from antennamod import antenna_pattern
import numpy as np
#import lalsimulation
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap


def sensitivity(t0,ra,dec,iota,psi,dist,m1inj,m2inj,det,fmaxhcode,massdet,ltotal,lsec):
    
    #load the data files:
    susp=np.loadtxt('datafiles/psus.txt')
    coat=np.loadtxt('datafiles/pcoat.txt')  
    seis=np.loadtxt('datafiles/pseismic.txt')
    quan=np.loadtxt('datafiles/pquantum.txt')
    i=0
    freqrange = coat[i*(500+1):(1+i)*500,0]   

#Calculate ISCO
    #masses in SI units:
    c = 299792458.        # m/s
    G = 6.67384e-11       # m**3 / s**2 / kg = M1+M2   # total mass
    M1 = m1inj * 1.9885469549614615e+30 #lal.MSUN_SI
    M2 = m2inj *  1.9885469549614615e+30#lal.MSUN_SI
    M = M1 + M2
    Mu = M1*M2/M    # reduced mass in CMS
    Mchirp = M**(2./5.)*Mu**(3./5.)   # chirp mass
    aorbI = 6.*G*M/c**2         # end inspiral at ISCO
    omI = np.sqrt(G*M/aorbI**3)
    forbI = omI/(2.*np.pi)
    fGWI = 2.*forbI
#Freqrange for h:
    ftheory = freqrange[freqrange < fmaxhcode] 
    ftheory = ftheory[ftheory < fGWI]    

    def htheory(m1inj,m2inj,dist,ftheory):
        #For optimal location and orientation:
        #masses in SI units:
        c = 299792458.        # m/s
        G = 6.67384e-11       # m**3 / s**2 / kg = M1+M2   # total mass
        M1 = m1inj * 1.9885469549614615e+30 #lal.MSUN_SI
        M2 = m2inj * 1.9885469549614615e+30 #lal.MSUN_SI
        M = M1 + M2
        Mu = M1*M2/M    # reduced mass in CMS
        dist = 1.e6* 3.085677581491367e+16 *dist #  lal.PC_SI*dist
        Mchirp = M**(2./5.)*Mu**(3./5.)   # chirp mass
        ht = np.sqrt(5*np.pi/24.)*G**2*Mchirp**2/(c**5*dist)* \
                np.power(np.pi*G*Mchirp*ftheory/c**3,-7./6.)
        return ht

    def tchirp(m1inj, m2inj, f):
        #For optimal location and orientation:
        #masses in SI units:
        c = 299792458.        # m/s
        G = 6.67384e-11       # m**3 / s**2 / kg = M1+M2   # total mass
        M1 = m1inj * 1.9885469549614615e+30 #lal.MSUN_SI
        M2 = m2inj * 1.9885469549614615e+30 #lal.MSUN_SI
        M = M1 + M2
        Mu = M1*M2/M    # reduced mass in CMS
        forb0 = f/2.
        fGW0 = f
        om0 = 2.*np.pi*forb0
        aorb0 = (G*M/om0**2)**(1./3.)       # Kepler's 3rd law
        vorb0 = aorb0*om0                   # orbital velocity in cms
        borb0 = vorb0/c                     # v/c
        #(2.19 en adV Detector)
        tchirp = 5./(256.*(Mu/M))*(G*M/c**3)/borb0**8   # chirp time
        return tchirp

    #Using antennaMatt insted 
#    def antenna_response( t0, ra, dec, psi, det ):
#        gps = lal.LIGOTimeGPS( t0 )
#        gmst_rad = lal.GreenwichMeanSiderealTime(gps)
#
#        # create detector-name map
#        detMap = {'H1': lal.LALDetectorIndexLHODIFF, \
#                    'H2': lal.LALDetectorIndexLHODIFF, \
#                    'L1': lal.LALDetectorIndexLLODIFF, \
#                    'G1': lal.LALDetectorIndexGEO600DIFF, \
#                    'V1': lal.LALDetectorIndexVIRGODIFF, \
#                    'T1': lal.LALDetectorIndexTAMA300DIFF}
#
#        try:
#            detector=detMap[det]
#        except KeyError:
#            raise ValueError, "ERROR. Key %s is not a valid detector name." % (det)
#
#          # get detector
#        detval = lal.CachedDetectors[detector]
#
#        response = detval.response
#
#        # actual computation of antenna factors
#        fp, fc = lal.ComputeDetAMResponse(response, ra, dec, psi, gmst_rad)
#        return fp, fc  
#
    ftotal =  antenna_pattern(det, np.radians(ra), np.radians(dec), psi, t0)

    if det == 0:
        det = 'H1'
    if det == 1:
        det = 'L1'
    if det == 2.:
        det = 'H1'
    
    #apt, act = antenna_response(t0, np.radians(ra), np.radians(dec), psi, det)
    htheory = htheory(m1inj,m2inj,dist,ftheory)
    
    # get GW strain at detector
    #h = htheory* np.sqrt(apt**2 + act**2)
    h = htheory * ftotal 

    #Get the noises:
    mm = np.arange(40.,220.,20.)
    
    #l4 = np.array([1, 2, 3])
    l4 = 100*np.array([0.6,0.85,1.1] )
    
    ltotallist = 100*np.array([1.6,1.87,2.14] )
    mm = int(np.where(mm == massdet)[0])
    ltotal = int(np.where(ltotallist == ltotal)[0])
    lsec = int(np.where(l4 == lsec)[0])
   

    if ltotal == 0:
        i=np.array([0,1,2])
    if ltotal == 1:
        i=np.array([3,4,5])
    if ltotal == 2:
        i=np.array([6,7,8])
    
    #print 'ltotal', ltotal
    mm= mm +1
    i = i[lsec]
    seis = seis[i*(500):(1+i)*500,mm]   
    coat = coat[i*(500):(1+i)*500,mm]   
    susp = susp[i*(500):(1+i)*500,mm]   
    quan = quan[i*(500):(1+i)*500,mm]

#ASD 
    asd = np.sqrt(seis**2 + coat**2 + susp**2 + quan**2)

##Extra ASD:
#    for i in np.arange(len(freqrange)):
#        asd[i] = lalsimulation.SimNoisePSDAdvVirgo(freqrange[i])
#        asd[i] = lalsimulation.SimNoisePSDaLIGONSNSOpt(freqrange[i])
#    asd = np.sqrt(asd)
#
    # get the SNR
    freqidx = freqrange[freqrange < 50]
    hidx =  h[0:len(freqidx)]
    idx = abs(abs(hidx*np.sqrt(freqrange[:len(hidx)])) - abs(asd[:len(hidx)]))
    idx = np.argmin(idx)

    idx2 = len(h)-1 
    if  fmaxhcode > 100 and fGWI > 100:
        freqidx2 = freqrange[freqrange > 60]
        idx2 = np.where(freqrange == freqidx2[0])[0][0]
        temp = abs(abs(h[idx2:]*np.sqrt(freqrange[idx2:len(h)])) - abs(asd[idx2:len(h)]))
        temp = np.argmin(temp)
        idx2 = temp + idx2

#For the Horizon Distance
    hwhite1 = htheory/asd[:len(htheory)] # whiten waveformi
    #deltaF = freqrange[2]-freqrange[1]
    deltaF =1.  #######PReguntar por esto
    snr1 = np.sqrt(4.*deltaF*np.vdot(hwhite1, hwhite1).real)

    #hwhite = h/asd[:len(h)] # whiten waveformi
    hwhite = h[idx:]/asd[idx:len(h)] # whiten waveformi
    hwhite = h[idx:idx2]/asd[idx:idx2] # whiten waveformi
    #deltaF = freqrange[2]-freqrange[1]
    deltaF = 1.
    snr = np.sqrt(4.*deltaF*np.vdot(hwhite, hwhite).real)

    timeelapsed =0.
    if snr > 2:
        timeelapsed = abs(tchirp(m1inj,m2inj,ftheory[idx])  - tchirp(m1inj,m2inj,ftheory[idx2]))
   
    #timeelapsed = abs(tchirp(m1inj,m2inj,ftheory[147])  - tchirp(m1inj,m2inj,ftheory[idx2]))
    
    time_to_ISCO =  tchirp(m1inj,m2inj,ftheory[-1])        
    BNS_horizon = snr1/8. *dist
    ave_range = BNS_horizon / 2.26
    #print 'i %d' %i
    #print 'mm %d' %mm

    return snr, h, asd, freqrange, ftheory, timeelapsed,seis, coat, \
            susp, fGWI, idx, idx2, quan, time_to_ISCO, BNS_horizon, ave_range
