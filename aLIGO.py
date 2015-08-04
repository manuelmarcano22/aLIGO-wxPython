#!/usr/bin/env python

"""
wxPython code to compute the SNR using LIGO data for Seismic, Coating and Suspension Noise.
Also computes the SNR for a given binary with equal masses at a particular distance. 
Y con calendario. 
"""

import numpy as np
import wx
from sensitivity import sensitivity
import os
import matplotlib
matplotlib.interactive(False)
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.figure import Figure
from matplotlib.pyplot import gcf, setp
from wx.lib.pubsub import pub
from mpl_toolkits.basemap import Basemap
class Knob:
    """
    Knob - simple class with a "setKnob" method.  
    A Knob instance is attached to a Param instance, e.g. param.attach(knob)
    Base class is for documentation purposes.
    """
    def setKnob(self, value):
        pass

class Param:
    """
    The idea of the "Param" class is that some parameter in the GUI may have
    several knobs that both control it and reflect the parameter's state, e.g.
    a slider, text, and dragging can all change the value of the frequency in
    the waveform of this example.  
    The class allows a cleaner way to update/"feedback" to the other knobs when 
    one is being changed.  Also, this class handles min/max constraints for all
    the knobs.
    """
    def __init__(self, initialValue=None, lista=np.linspace(0,1,1000)):
        self.minimum = round(min(lista))
        self.maximum = round(max(lista))
        self.lista = lista
        self.value = initialValue
        self.knobs = []
        
    def attach(self, knob):
        self.knobs += [knob]
        
    def set(self, value, knob=None):
        self.value = value
        self.value = self.constrain(value)
        for feedbackKnob in self.knobs:
            if feedbackKnob != knob:
                feedbackKnob.setKnob(self.value)
        return self.value


    def constrain(self, value):
        cercano = min(self.lista, key=lambda x:abs(x-value))
        value = cercano
        return value

class Param2:
#    Param pero Sin lo de Constrain
    def __init__(self, initialValue=None, lista=np.linspace(0,1,1000)):
        self.minimum = round(min(lista))
        self.maximum = round(max(lista))
        self.lista = lista
        self.value = initialValue
        self.knobs = []
        
    def attach(self, knob):
        self.knobs += [knob]
    
    def set(self, value, knob=None):
        self.value = value
        for feedbackKnob in self.knobs:
            if feedbackKnob != knob:
                feedbackKnob.setKnob(self.value)
        return self.value
        
class SpinerGroup(Knob):
    def __init__(self, parent, position, title, label,size,param):
        self.spinertitle = wx.StaticText(parent, label=title)
        self.spinerlabel = wx.StaticText(parent, label=label)


        pub.subscribe(self._update, 'valuechanged') 
            
        self.spiner = wx.SpinCtrlDouble(parent,-1,pos=wx.DefaultPosition,size=size)
        
        self.spiner.SetRange(min(param.lista),max(param.lista))
        self.setKnob(param.value)
        self.spiner.Bind(wx.EVT_SPINCTRLDOUBLE,self.spinerHandler)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.spinertitle,0)
        sizer.Add(self.spiner, 1, wx.EXPAND)
        self.sizer = sizer
        
        self.param = param
        self.param.attach(self)
        self.fourierDemoWindow = FourierDemoWindow
    
    def spinerHandler(self, evt):
        value = evt.GetValue() #/ 10000.
        self.param.set(value)

    def setKnob(self, value):
        self.spiner.SetValue(value)

    def _update(self, maxfreq):
        self.spiner.SetMax(int(maxfreq))

class SliderGroup(Knob):
    def __init__(self, parent, label, param,orientation):
        self.sliderLabel = wx.StaticText(parent, label=label)
        
        font = wx.Font(12, wx.MODERN, wx.ITALIC, wx.NORMAL)
        
        self.sliderText = wx.TextCtrl(parent, -1, style=wx.TE_PROCESS_ENTER)
        self.sliderLabel.SetFont(font)

        self.slider = wx.Slider(parent, -1,style=orientation   )
        self.slider.SetMax(max(param.lista))
        self.slider.SetMin(min(param.lista))
        self.setKnob(param.value)
        
        sizer = wx.BoxSizer(orientation)
        
        sizer.Add(self.sliderLabel, 0, wx.EXPAND | wx.ALIGN_CENTER | wx.ALL, border=2)
        sizer.Add(self.sliderText, 0, wx.EXPAND | wx.ALIGN_CENTER | wx.ALL, border=2)
        sizer.Add(self.slider, 1, wx.EXPAND)
        self.sizer = sizer

        self.slider.Bind(wx.EVT_SLIDER, self.sliderHandler)

        self.sliderText.Bind(wx.EVT_TEXT_ENTER, self.sliderTextHandler)

        self.param = param
        self.param.attach(self)

    def sliderHandler(self, evt):
        value = evt.GetInt() #/ 10000.
        self.param.set(value)
        
    def sliderTextHandler(self, evt):
        value = float(self.sliderText.GetValue())
        self.param.set(value)
        
        
    def setKnob(self, value):
        self.sliderText.SetValue('%g'%value)
        #self.slider.SetValue(value*10000)
        self.slider.SetValue(value)


class FourierDemoFrame(wx.Frame):
    def __init__(self, *args, **kwargs):
        wx.Frame.__init__(self, *args, **kwargs)

        self.create_menu()
        self.statusbar=self.CreateStatusBar(style=wx.SB_RAISED)
        self.statusbar.SetFieldsCount(2)
        self.statusbar.SetStatusWidths([-3, -3])
        pub.subscribe(self.change_statusbar, 'changestatusbar')

        self.fourierDemoWindow = FourierDemoWindow(self)
        self.frequencySliderGroup = SliderGroup(self, label='mass det (kg):', \
            param=self.fourierDemoWindow.f0,orientation =  wx.SL_HORIZONTAL )
        self.amplitudeSliderGroup = SliderGroup(self, label='L4 (cm)', \
                param=self.fourierDemoWindow.A, orientation = wx.SL_HORIZONTAL)
        
        self.ltotalSliderGroup = SliderGroup(self, label='Total Length Susp (cm)', \
                param=self.fourierDemoWindow.ltotal, orientation = wx.SL_HORIZONTAL)
        
        self.massSliderGroup = SliderGroup(self, label='Solar Mass', \
            param=self.fourierDemoWindow.mass,orientation = wx.SL_VERTICAL)

        self.distSliderGroup = SliderGroup(self, label='Dist (Mpc)', \
            param=self.fourierDemoWindow.dist,orientation = wx.SL_VERTICAL)
    
        self.radioBoxes = RadioBoxFram(self, 'Detectors',['H1','L1'], \
                detector = self.fourierDemoWindow.detector)
        
        self.spinerFreq = SpinerGroup(self,title='Maximum Waveform Freq (Hz)',label='l',\
                position=wx.DefaultPosition, size=wx.DefaultSize,\
                param=self.fourierDemoWindow.freqmax)

        sizerbox = wx.BoxSizer(wx.HORIZONTAL)
        sizerbox.Add(self.massSliderGroup.sizer, 0  ,\
                wx.EXPAND|wx.ALIGN_CENTER | wx.ALL, border=5)
        
        sizerbox.Add(self.distSliderGroup.sizer, 0 ,\
                wx.EXPAND | wx.ALIGN_CENTER | wx.ALL, border=5)

        sizer3 = wx.BoxSizer(wx.VERTICAL)
        sizer3.Add(sizerbox,1,wx.EXPAND)
        sizer3.Add(self.spinerFreq.sizer,0)
        sizer3.Add(self.radioBoxes.radioBox,0,wx.EXPAND)

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.fourierDemoWindow,1, wx.EXPAND)
        sizer.Add(sizer3,0,wx.EXPAND)

        sizer2 = wx.BoxSizer(wx.VERTICAL)
        sizer2.Add(sizer,1)
        
        sizer2.Add(self.frequencySliderGroup.sizer, 0, \
                wx.EXPAND | wx.ALIGN_CENTER | wx.ALL, border=5)
        sizer2.Add(self.amplitudeSliderGroup.sizer, 0, \
                wx.EXPAND | wx.ALIGN_CENTER | wx.ALL, border=5)
        sizer2.Add(self.ltotalSliderGroup.sizer,0,\
                wx.EXPAND | wx.ALIGN_CENTER | wx.ALL, border=5)

        self.SetSizer(sizer2)
        

    def change_statusbar(self, msg1,msg2):
        pass
        self.statusbar.SetStatusText(msg1,0)
        self.statusbar.SetStatusText(msg2,1)



    def create_menu(self):
        self.menubar = wx.MenuBar()
        
        menu_file = wx.Menu()
        m_expt = menu_file.Append(-1, "&Save plot\tCtrl-S", "Save plot to file")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_expt)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)
        
        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&About\tF1", "About the demo")
        self.Bind(wx.EVT_MENU, self.on_about, m_about)
       
        menu_plot = wx.Menu()
        
        m_map = menu_plot.AppendRadioItem(-1, 'Blue Marble')
        m_antenna = menu_plot.AppendRadioItem(-1, 'Antenna Pattern')
        self.Bind(wx.EVT_MENU, self.on_antenna, m_antenna)
        self.Bind(wx.EVT_MENU, self.on_map, m_map)

        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_plot, "&Plot")
        self.menubar.Append(menu_help, "&Help")
        self.SetMenuBar(self.menubar)


    def on_save_plot(self, event):
        file_choices = "PNG (*.png)|*.png"
        
        dlg = wx.FileDialog(
            self, 
            message="Save plot as...",
            defaultDir=os.getcwd(),
            defaultFile="plot.png",
            wildcard=file_choices,
            style=wx.SAVE)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.fourierDemoWindow.canvas.print_figure(path, dpi=10000)
            #self.flash_status_message("Saved to %s" % path)
        
    def on_exit(self, event):
        self.Destroy()
        
    def on_about(self, event):
        msg = """ A demo using wxPython with matplotlib:
         * Calculate the Sensitivity Curve iLIGO
         * Save the plot to a file using the File menu
         * It is Stationary phased approximation
         * L4 is the seccion Lenght
         * LTotal
        """
        dlg = wx.MessageDialog(self, msg, "About", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()

    def on_antenna(self, event):
        pub.sendMessage('detmapchanged' , newmap =1. )#,newdet = 0. ,newmap = 1.)      
    
    def on_map(self, event):
        pub.sendMessage('detmapchanged' , newmap = 0.)#,newdet =0.  ,newmap = 0.)      


class RadioBoxFram(Knob):
    def __init__(self, parent, label, paramlist, detector ):
        
        self.radioBox = wx.RadioBox(parent, -1, label, (10,10), wx.DefaultSize, \
                paramlist, 2, wx.RA_SPECIFY_COLS)
    
        self.detector = detector 
        self.detector.attach(self) 
        
        self.radioBox.Bind(wx.EVT_RADIOBOX, self.radioHandler)

    def radioHandler(self, evt):
        value = self.radioBox.GetSelection()
        value = float(value)
        self.detector.set(value)
        pub.sendMessage('detmapchanged', newdet = 1.)      

class FourierDemoWindow(wx.Window, Knob):
    def __init__(self, *args, **kwargs):
        wx.Window.__init__(self, *args, **kwargs)
        self.lines = []
        self.figure = Figure()
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.state = ''
        
        
        self.aaprueba =   Param2(0., lista = np.array([0.,1.]))
        self.aaprueba.attach(self)
        
       
        self.mapastatus =   Param(0., lista = np.array([0.,1.]))   

       #mass det
        self.f0 =   Param(40., lista=np.arange(40,220,20 ) )
        
        #L4 
        self.A = Param(0.85*100, lista = 100*np.array([0.6,0.85,1.1] ))
        
        #ltotal
        self.ltotal = Param(1.6*100, lista = 100*np.array([1.6,1.87,2.14]))

        #FOr Mass and Distance
        self.mass = Param(4., lista =np.arange(1,11,.1) )
        self.dist = Param(100.,  lista=np.arange(10,10000,10)) 
       

        self.detector = Param(0., lista=np.array([0.,1.]) )
        
        self.dec = Param(-29.4,lista=np.linspace(-90 ,90,180))
        self.ra = Param(71.85-180, lista=np.linspace(-180,180,360) )

        self.ra.attach(self)
        self.dec.attach(self)
        
        #Para detector
        self.f0.attach(self)
        self.A.attach(self)
        self.mass.attach(self)
        self.dist.attach(self)
        self.detector.attach(self)
        self.ra.attach(self)
        self.dec.attach(self)
        self.mapastatus.attach(self)
        self.ltotal.attach(self)
        self.Bind(wx.EVT_SIZE, self.sizeHandler)

        self.freqmax = Param(100, lista= np.arange(11,2000,1) ) 
        self.freqmax.attach(self)
        self.draw()
        pub.subscribe(self._update, 'detmapchanged') 
    def _update(self, newdet=None, newmap = None, fechac =None, horac = None   ):
            
        mapita = Basemap(ax=self.inset_ax,projection='moll',lon_0 = 0)

        X = np.arange(-180,180.,1.)
        Y = np.arange(-90.,90.,1.)
        X, Y = np.meshgrid(X,Y)
                
        if newmap == 1.:
            #mapita.ax.clear() 
            
            self.mapastatus.set(1.)
            self.aaprueba.set(50.)
            
            if self.detector.value == 0:
                Z = np.loadtxt('datafiles/H1antenna.out')
            
            if self.detector.value == 1:
                Z = np.loadtxt('datafiles/L1antenna.out')
           
            im1 = mapita.pcolormesh(X,Y,Z, shading='flat',\
                    cmap=matplotlib.cm.jet,latlon=True)
            self.scatter = mapita.scatter(self.raplot, self.decplot,c='r',marker='*',s=100)
            mapita.drawcoastlines()

        if newmap ==0.:
            self.mapastatus.set(0.)
            self.aaprueba.set(10.)
            self.scatter = mapita.scatter(self.raplot, self.decplot,c='r',marker='*',s=100)
            mapita.ax.clear() 
            mapita.bluemarble(scale=.1)
            self.scatter = mapita.scatter(self.raplot, self.decplot,c='r',marker='*',s=100)

        
        if newdet == 1.:

            if self.mapastatus.value == 1.:
                if self.detector.value == 0:
                    Z = np.loadtxt('datafiles/H1antenna.out')
                
                if self.detector.value == 1:
                    Z = np.loadtxt('datafiles/L1antenna.out')
               
                #mapita.ax.clear() 
                im1 = mapita.pcolormesh(X,Y,Z, shading='flat',\
                        cmap=matplotlib.cm.jet,latlon=True)
                self.scatter = mapita.scatter(self.raplot, self.decplot,c='r',marker='*',s=100)
                mapita.drawcoastlines()
            
        self.repaint()

    def sizeHandler(self, *args, **kwargs):
        self.canvas.SetSize(self.GetSize())
        

    def OnLeftDown(self, event):
        #pt = event.xdata  # position tuple
        mapita3 = Basemap(projection='moll',lon_0 = 10,resolution=None)
        valuex = (event.xdata)#*(180./np.pi)
        valuey = (event.ydata)#*(180./np.pi)
        valuex, valuey = mapita3(valuex,valuey ,inverse=True )
        self.ra.set(valuex)
        self.dec.set(valuey)
        #print 'ra and dec',self.ra.value +180. , self.dec.value

    def draw(self):
        if not hasattr(self, 'subplot1'):
            self.subplot1 = self.figure.add_subplot(111)
            self.inset_ax = self.figure.add_axes([0.55,0.04,0.45,0.3] )
            mapita = Basemap(ax=self.inset_ax,projection='moll',lon_0 = 0,resolution=None)

            # hook some mouse events
            self.canvas.callbacks.connect('button_press_event', self.OnLeftDown)
        
        #Call to Functionn
        t0=900000000 #gps seconds
        iota = 0. #(radians)
        psi = 0.343 #(radians)
        snr1, h1, asd1, freqrange1, freqrangetheory, timee1,\
                seis1, coat1, susp1, fGWI, idx1,idx2, quantum1 = \
                self.computenoise(t0, self.ra.value+180., \
                self.dec.value, iota, psi, self.dist.value, \
                self.mass.value, self.mass.value, \
                self.detector.value, self.freqmax.value, \
                self.f0.value, self.ltotal.value, self.A.value)

        
        self.raplot, self.decplot = mapita(self.ra.value, self.dec.value)
        mapita.bluemarble(scale=.2)
#Map
        ###
        self.subplot1.set_yscale('log')
        self.subplot1.set_xscale('log')

        #Text box for SNR
        annotation_string = r'$SNR=%.2f$'%snr1
        annotation_string +='\n'
        annotation_string +=r'$t_{elapsed}=%.2f$'%timee1
        self.text = self.subplot1.text(0.02,0.97, annotation_string,size=20 ,ha='left', va='top',\
                transform=self.subplot1.transAxes, \
                bbox = dict(boxstyle='round', ec=(1., 0.5, 0.5),\
                            fc=(1., 0.8, 0.8),))

        self.subplot1.set_ylim([1e-27, 1e-17])
        self.subplot1.set_xlim([1, 2e4])

        self.subplot1.grid(True, which="both", ls="-", alpha=.5)

        self.subplot1.set_xlabel(r"$frequency, f(Hz)$", fontsize=18)
        self.subplot1.set_ylabel(r"$S^{1/2}_h (f) \, (Hz^{-1/2})$",fontsize=18)

#coating Brownian, seismic, suspension thermal noise.

        self.lines += self.subplot1.plot(freqrange1, asd1, \
                '-', linewidth = 5,label='Sensitivity Curve')
        self.lines += self.subplot1.plot(freqrange1,seis1,\
                '-',label='Seismic Noise')
        self.lines += self.subplot1.plot(freqrange1,susp1,\
                '-', label='Suspension Thermal Noise')
        self.lines += self.subplot1.plot(freqrange1,coat1,\
                '-', label='Coating Brownian Noise')
        
        
        self.lines += self.subplot1.plot(freqrangetheory,h1*np.sqrt(freqrangetheory),\
                '-', label='Effective Induced Strain')
        
        self.lines += self.subplot1.plot(freqrange1,np.sqrt(quantum1),\
                '-', label='Quantum Noise')
        
        self.scatter = mapita.scatter(self.raplot, self.decplot,c='r', marker = '*', s=100)
       
        if snr1 > 0.5:
            self.scatteridx = self.subplot1.scatter(freqrangetheory[idx1],\
                    h1[idx1]*np.sqrt(freqrangetheory[idx1]),c='r',s=100)
            #self.scatteridx2 = self.subplot1.scatter(freqrangetheory[idx2],\
            #            h1[idx2]*np.sqrt(freqrangetheory[idx2]),c='r',s=100)
        #print idx1 

        self.subplot1.legend(loc='best', fancybox=True, framealpha=0.5)
      
        mapita.bluemarble(scale = .1)

#Set some plot attributes

        self.subplot1.set_title("Click and drag sliders to change source and design parameters", fontsize=14)

    def computenoise(self,t0,ra,dec,iota,psi,dist,m1inj,m2inj,det,fmaxhcode,massdet,ltotal,lsec):
        a, b, c, d, e,f,g,h,i,j,k,l,m = sensitivity(t0,ra,dec,iota,psi,dist,m1inj,m2inj,det,fmaxhcode,massdet,ltotal,lsec)
        return a, b, c, d, e, f,g,h,i, j,k,l,m
    
    def repaint(self):
        self.canvas.draw()

    def setKnob(self, value):
        # Note, we ignore value arg here and just go by state of the params
       
       #Call to Function
        t0=900000000 #gps seconds
        iota = 0. #(radians)
        psi = 0.343 #(radians)

        snr1, h1, asd1, freqrange1, freqrangetheory, timee1, \
                seis1, coat1, susp1, fGWI,idx1,idx2, quantum1 = \
                self.computenoise(t0, self.ra.value+180., self.dec.value, iota, psi, self.dist.value, \
                self.mass.value, self.mass.value, self.detector.value, self.freqmax.value, \
                self.f0.value, self.ltotal.value, self.A.value)
        pub.sendMessage('valuechanged',maxfreq = fGWI) # freqrange1[len(timesall1)])      

        annotation_string = r'$SNR=%.2f$'%snr1
        annotation_string +='\n'
        annotation_string +=r'$t_{elapsed}=%.2f$'%timee1
        self.text.set_text(annotation_string  )

        #prueba = 125 #nip.where(freqrange1== 10.046251)
        #print 'freq', freqrange1[prueba]
        #print 'seismic:,', seis1[prueba]
        #print 'coat:', coat1[prueba]
        #print 'Susp::', susp1[prueba]

        setp(self.lines[0], xdata=freqrange1, ydata=asd1)
        setp(self.lines[1], xdata=freqrange1, ydata=seis1)
        setp(self.lines[2], xdata=freqrange1, ydata=susp1)
        setp(self.lines[3], xdata=freqrange1, ydata=coat1)
        setp(self.lines[4], xdata=freqrangetheory, ydata=h1*np.sqrt(freqrangetheory))
        setp(self.lines[5], xdata=freqrange1, ydata=np.sqrt(quantum1))
        
        mapita2 = Basemap(projection='moll',lon_0 = 0,resolution=None)
        self.raplot, self.decplot = mapita2(self.ra.value, self.dec.value)
        self.scatter.set_offsets([ self.raplot, self.decplot])
        
        if snr1 < 0.5:
            self.scatteridx.set_offsets([0,0])
            #self.scatteridx2.set_offsets([0,0])
        if snr1 > 0.5:
            self.scatteridx.set_offsets([freqrangetheory[idx1],\
                    h1[idx1]*np.sqrt(freqrangetheory[idx1])])
#            self.scatteridx2.set_offsets([freqrangetheory[idx2],\
#                    h1[idx2]*np.sqrt(freqrangetheory[idx2])])
#            if max(freqrangetheory) < 200:
#                self.scatteridx2.set_offsets([0,0])

        #print self.raplot
        if self.detector.value == 0:
            sbdet = 'Harford'
        if self.detector.value ==1:
            sbdet ='Livinston'

        pub.sendMessage('changestatusbar',\
                msg1= \
                'Horizon Distance for SNR of 8:    %s'% 'algo' +'    Mpc' ,\
                msg2 = 'Time to ISCO: %s' %'algo' + ' s')
#                msg2='Detector: %s' % sbdet)
        
        self.repaint()

class App(wx.App):
    def OnInit(self):
        self.frame1 = FourierDemoFrame(parent=None, title="aLIGO Demo", size=(900, 900))  #640, 480
        self.frame1.Show()
        return True

app = App()
app.MainLoop()
