#!/usr/bin/env python
import lal
import lalsimulation
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#from scipy import interpolate


def sensitivity(t0,ra,dec,iota,psi,dist,m1inj,m2inj,det,fmaxhcode,massdet,ltotal,lsec):
    
    #load the data files:
    susp=np.loadtxt('datafiles/psus.txt')
    coat=np.loadtxt('datafiles/pcoat.txt')  
    seis=np.loadtxt('datafiles/pseismic.txt')
    i=0
    freqrange = coat[i*(500+1):(1+i)*500,0]   

#Calculate ISCO
    #masses in SI units:
    c = 299792458.        # m/s
    G = 6.67384e-11       # m**3 / s**2 / kg = M1+M2   # total mass
    M1 = m1inj * lal.MSUN_SI
    M2 = m2inj * lal.MSUN_SI
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
        M1 = m1inj * lal.MSUN_SI
        M2 = m2inj * lal.MSUN_SI
        M = M1 + M2
        Mu = M1*M2/M    # reduced mass in CMS
        dist = 1.e6*lal.PC_SI*dist
        Mchirp = M**(2./5.)*Mu**(3./5.)   # chirp mass
        ht = np.sqrt(5*np.pi/24.)*G**2*Mchirp**2/(c**5*dist)* \
                np.power(np.pi*G*Mchirp*ftheory/c**3,-7./6.)
        return ht

    def tchirp(m1inj, m2inj, f):
        #For optimal location and orientation:
        #masses in SI units:
        c = 299792458.        # m/s
        G = 6.67384e-11       # m**3 / s**2 / kg = M1+M2   # total mass
        M1 = m1inj * lal.MSUN_SI
        M2 = m2inj * lal.MSUN_SI
        M = M1 + M2
        Mu = M1*M2/M    # reduced mass in CMS
        forb0 = f/2.
        fGW0 = f
        om0 = 2.*np.pi*forb0
        aorb0 = (G*M/om0**2)**(1./3.)       # Kepler's 3rd law
        vorb0 = aorb0*om0                   # orbital velocity in cms
        borb0 = vorb0/c                     # v/c
        tchirp = 5./(256.*(Mu/M))*(G*M/c**3)/borb0**8   # chirp time
        return tchirp
    
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
        return fp, fc  

    if det == 0:
        det = 'H1'
    if det == 1:
        det = 'L1'
    
    apt, act = antenna_response(t0, np.radians(ra), np.radians(dec), psi, det)
    htheory = htheory(m1inj,m2inj,dist,ftheory)
    # get GW strain at detector
    h = htheory* np.sqrt(apt**2 + act**2)

    #Get the noises:
    mm = np.arange(40.,220.,20.)
    
    #l4 = np.array([1, 2, 3])
    l4 = 100*np.array([0.6,0.85,1.1] )
    
    #ltotallist = np.array([1,2,3])
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
    i = i[lsec]
    mm= mm +1
    seis = seis[i*(500):(1+i)*500,mm]   
    coat = coat[i*(500):(1+i)*500,mm]   
    susp = susp[i*(500):(1+i)*500,mm]   
    #asd = np.sqrt(seis**2 + coat**2 + susp**2)
    asd = np.sqrt(seis**2 + coat**2 + susp**2)

    
    #Interpolation:
    #hinter = interpolate.interp1d(ftheory, h)
    #asdinter = interpolate.interp1d(ftheory, asd[:len(h)])
    #xnew = np.linspace(min(ftheory), max(ftheory), 50000)
    #hinter = hinter(xnew)
    #asdinter = asdinter(xnew)

    # get the SNR
    #idx = np.argwhere(np.isclose(hinter,asdinter[len(h)])).reshape(-1) 
    freqidx = freqrange[freqrange < 100]
    hidx =  h[0:len(freqidx)]
    idx = abs(abs(hidx*np.sqrt(freqrange[:len(hidx)])) - abs(asd[:len(hidx)]))
    #idx = np.where(np.sort(idx)==idx)
    #idx = idx[0:50]
#    idx = np.sort(idx)
#    idx = np.where(idx == idx[-1])[0][0]
    idx = np.argmin(idx)

##    idxtemp=np.array([])
#    for i in range(len(h)-1):
#        xnew = np.linspace(ftheory[i],ftheory[i+1],1000)
#        yinterph = np.array(np.interp(xnew, ftheory[i:i+1], h[i:i+1]))
#        yinterpasd = np.array(np.interp(xnew, ftheory[i:i+1], asd[i:i+1]))
#        idx2 = np.where(yinterpasd == yinterph)
#        if np.sum(idx2) > 0:
#            np.append(idxtemp,i)
#    print idxtemp
#
#
##    
#    idx2 = abs(abs(h*np.sqrt(freqrange[:len(h)])) - abs(asd[:len(h)]))
#    idx2=zip(*[freqrange,idx2])
#    idx2desde = idx2[idx2[:,0]>1000]
#    idx2desde = np.sort(idx2desde, axis=1) 
#    
#    
    idx2 = 0
    desde = 200.
    if max(ftheory) > desde:
        freqidx2 = ftheory[ftheory > desde]
        desdeidx2 = np.where( ftheory == freqidx2[0])[0][0]
        hastaidx2 = np.where(ftheory == freqidx2[-1])[0][0]
        #print 'freqdedsde', ftheory[desdeidx2], ftheory[hastaidx2]
        hidx2 =  h[desdeidx2:hastaidx2]
        idx2 = abs(abs(hidx2*np.sqrt(ftheory[desdeidx2:hastaidx2])) - abs(asd[desdeidx2:hastaidx2]))
        idx2 = np.argmin(idx2)
        #print 'idx2, freq, idx1', idx2, freqrange[idx2], idx
    #prine idx
    #print 'h , asd', freqrange[idx],h[idx],asd[idx] 
    #print idx[0:10]
    #cercano = min(h, key=lambda x:abs(x-value))
    #idx = np.where(hinter==asdinter)
    #i=  idx#np.where(freqrange > 15 )
    #print i
    #idx=1
    hwhite = h[idx:]/asd[idx:len(h)] # whiten waveformi
    deltaF = freqrange[2]-freqrange[1]
    #print deltaF
    snr = np.sqrt(4.*deltaF*np.vdot(hwhite, hwhite).real)

    timeelapsed = abs(tchirp(m1inj,m2inj,ftheory[0])  - tchirp(m1inj,m2inj,ftheory[-1]))

    #print 'i %d' %i
    #print 'mm %d' %mm

    return snr, h, asd, freqrange, ftheory, timeelapsed, seis, coat, susp, fGWI, idx, idx2
#
##
##Prueba:s
##To call the FUnctions:
#t0=900000000 #gps seconds
#iota = 0. #(radians) Inclination
#psi = 0.343 #(radians) Polarization angle
#ra = 35.
#dec = 35.
#dist = 300.
#m1inj = 1
#m2inj = 1
#det =0
#fmaxhcode = 5000000.
#massdet = 80.
#ltotal =  100*2.14 # 1.6#1.87# 2.14 #1.87   # 2.1
#lsec = 100*1.1 #0.6
#
#snr1, h1, asd1, freqrange1, freqrangetheory, timee1, seis1, coat1, susp1, fGWI ,idx,idx2= \
#        sensitivity(t0,ra,dec,iota,psi,dist,m1inj,m2inj,det,fmaxhcode,massdet,ltotal,lsec)
#
#figure = plt.figure()
###
#subplot1 = figure.add_subplot(111)
#
#inset_ax = figure.add_axes([0.55,0.04,0.45,0.3] )
#mapita = Basemap(projection='moll',lon_0 = 10)
##, projection =ccrs.PlateCarree())
####
#subplot1.set_yscale('log')
#subplot1.set_xscale('log')
#
##Text box for SNR
#annotation_string = r'$SNR=%d$'%snr1
#annotation_string +='\n'
#annotation_string +=r'$t_{elapsed}=%f$'%timee1
#text = subplot1.text(0.02,0.97, annotation_string,size=20 ,ha='left', va='top',\
#        transform=subplot1.transAxes, \
#        bbox = dict(boxstyle='round', ec=(1., 0.5, 0.5),\
#                    fc=(1., 0.8, 0.8),))
#
#subplot1.set_ylim([1e-26, 1e-10])
#subplot1.set_xlim([1, 2e4])
#
#subplot1.grid(True, which="both", ls="-", alpha=.5)
#
#subplot1.set_xlabel(r"$frequency, f(Hz)$", fontsize=18)
#subplot1.set_ylabel(r"$S^{1/2}_h (f) \, (Hz^{-1/2})$",fontsize=18)
#
##coating Brownian, seismic, suspension thermal noise.
#
#
#lines=[]
#
#
#freqrange = np.arange(1,10001,1)
#q=np.zeros(len(freqrange)+1)
#
#for i in freqrange:
#        temp = lalsimulation.SimNoisePSDaLIGOQuantumZeroDetHighPower(i)
#        q[i] = temp
#
#subplot1.plot(freqrange, np.sqrt(q[:len(freqrange)]))
#
#
#subplot1.plot(freqrange1, asd1, \
#        '-', label='Sensitivity Curve')
#subplot1.plot(freqrange1,seis1,\
#        '-',label='Seismic Noise')
#subplot1.plot(freqrange1,susp1,\
#        '-', label='Suspension Thermal Noise')
#subplot1.plot(freqrange1,coat1,\
#        '-', label='Coating Brownian Noise')
#
#subplot1.plot(freqrangetheory,h1*np.sqrt(freqrangetheory),\
#        '-', label='Effective Induced Strain')
#
#subplot1.legend(loc='best', fancybox=True, framealpha=0.5)
#
#try:
#    subplot1.scatter(freqrangetheory[idx],h1[idx]*np.sqrt(freqrangetheory[idx]),c='r',s=500,marker='*')
#    subplot1.scatter(freqrangetheory[idx2],h1[idx2]*np.sqrt(freqrangetheory[idx2]),color='b',s=500)
#except:
#    print 'no se pudo'
##mapita.drawcoastlines()
#
#x,y = mapita(ra,dec)
#mapita.scatter(x,y,color='r')
#mapita.bluemarble(scale=.2)
#inset_ax = figure.add_axes([0.55,0.04,0.45,0.3] )
#
##Antenna Pathern:
#def antenna_response( t0, ra, dec, psi, det ):
#    gps = lal.LIGOTimeGPS( t0 )
#    gmst_rad = lal.GreenwichMeanSiderealTime(gps)
#
#    # create detector-name map
#    detMap = {'H1': lal.LALDetectorIndexLHODIFF, \
#    'H2': lal.LALDetectorIndexLHODIFF, \
#    'L1': lal.LALDetectorIndexLLODIFF, \
#    'G1': lal.LALDetectorIndexGEO600DIFF, \
#    'V1': lal.LALDetectorIndexVIRGODIFF, \
#    'T1': lal.LALDetectorIndexTAMA300DIFF}
#
#    try:
#        detector=detMap[det]
#    except KeyError:
#        raise ValueError, "ERROR. Key %s is not a valid detector name." % (det)
#
#    # get detector
#    detval = lal.CachedDetectors[detector]
#
#    response = detval.response
#
#    # actual computation of antenna factors
#    fp, fc = lal.ComputeDetAMResponse(response, ra, dec, psi, gmst_rad)
#    ft = np.sqrt(fp**2 + fc**2)
#    return ft
#t0=900000000 #gps seconds
#psi = 0#0.343 #(radians) Polarization angle
#det = 'L1'
#
#
#Xlista = np.arange(0,360.,1.)
#Ylista = np.arange(-90.,90.,1.)
#Z = np.array([])
#
#
#for i in Ylista:
#    for j in Xlista:
#        ztemp = antenna_response(t0, np.radians(j),np.radians(i), psi, det ) 
#        Z = np.append(Z,ztemp)
#
#X = np.arange(-180,180.,1.)
#Y = np.arange(-90.,90.,1.)
#X, Y = np.meshgrid(X,Y)
#Z = np.reshape(Z, X.shape)
###np.savetxt('L1antenna.out',Z)
###Z = np.loadtxt('datafiles/L1antenna.out')
###Contour
#im1 = mapita.pcolormesh(X,Y,Z,shading='flat',cmap=plt.cm.jet,latlon=True)
#mapita.drawcoastlines()
##
###For back to bluemarble
###xnot = np.zeros([5,5])
###im1.set_array(xnot)
###mapita.bluemarble(scale=.2)
##
##
##
##
###, projection =ccrs.PlateCarree())
###scatter = inset_ax.scatter(ra, dec,c='r',transform=ccrs.PlateCarree())
###scatter = inset_ax.scatter(0,0,transform=ccrs.PlateCarree()) 
###inset_ax.coastlines()
###inset_ax.gridlines(crs=ccrs.PlateCarree())
###        ,transform=ccrs.PlateCarree())
###inset_ax.plot(ra,dec,transform=ccrs.PlateCarree())
###fre10 = int(np.where(freqrange1 == 10.)[0])
###print 'at %d  is %d' %fre10 %susp1[fre10]
###figure = plt.figure()
###subplot2 = figure.add_subplot(111)
###subplot2.plot(freqrangetheory,h1*np.sqrt(freqrangetheory),\
###        '-', label='Effective Induced Strain')
###subplot2.set_yscale('log')
###subplot2.set_xscale('log')
###        
###        True , which="both", ls="-", alpha=.5)
#plt.show()
###

