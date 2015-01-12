''' graph with parameters that can 
be changed using sliders
Jan 2015'''

from numpy import (sin,arange,pi,sqrt,exp,arctan)
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

def homfun(t):
    '''homogeneous solution
    returns a function'''
    w=2*pi*sfreq.val
    tau=stau.val
    return 1/tau/sqrt(w**2+1/tau**2)*\
    w/sqrt(w**2+1/tau**2)*exp(-t/tau)
   
def partfun(t):
    '''particular solution
    returns a function'''
    w=2*pi*sfreq.val
    tau=stau.val
    return 1/tau/sqrt(w**2+1/tau**2)*\
    sin(w*t-arctan(w*tau))

###Plotting
plt.figure(1)
plt.clf()

#create sliders
axtau = plt.axes([.1, .025, 0.8, 0.025])
axfreq  = plt.axes([.1, .075, 0.8, 0.025])
#the values after the label are 
#the min, max, and initial values
stau = Slider(axtau, 'tau', 0.01, 2.0, valinit=1)
sfreq = Slider(axfreq, 'freq', 0.01, 2.0, valinit=1)
    
#create initial plot   
axplt = plt.axes([.1, .15, .8, .8])
tmax=7.
t = arange(0.0, tmax, tmax/500)
s = sin(2*pi*sfreq.val*t)
inp, = plt.plot(t,s, lw=1, color='green')
hom, = plt.plot(t,homfun(t), lw=2, color='red')
part, = plt.plot(t,partfun(t), lw=2, color='blue')
tot, = plt.plot(t,homfun(t)+partfun(t), lw=3, color='black')
plt.axis([0,tmax,-1.1,1.1])
plt.grid(which='both')

def update(val):
    '''draws the lines based on slider values'''
    w = 2*pi*sfreq.val
    #tau=stau.val
    inp.set_ydata(sin(w*t))
    hom.set_ydata(homfun(t))
    part.set_ydata(partfun(t))
    tot.set_ydata(homfun(t)+partfun(t))
    #fig.canvas.draw_idle()
sfreq.on_changed(update)
stau.on_changed(update)

#plt.show()
