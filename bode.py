'''basic bode plot
Jan 2015'''

from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
pi=np.pi


#create the transfer function
h = signal.lti([1], [1,1])
#set the range
w=np.logspace(-2,2)
#compute the magnitude and phase
w, mag, phase = signal.bode(h,w)


###Plotting
plt.figure(1)
plt.clf()

plt.subplot(2,1,1)
plt.loglog(w, 10**(mag/20),lw=2)
plt.grid(which='both')
plt.ylabel('magnitude')
plt.margins(.1)

plt.subplot(2,1,2)
plt.semilogx(w, phase/180,lw=2)
plt.grid(which='both')
plt.ylabel('phase/pi')
plt.xlabel('$w$ (rad/s)')
plt.margins(.1)

plt.show()

