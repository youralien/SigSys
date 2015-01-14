from numpy.fft import (fft,fftshift)
from matplotlib.pyplot import (plot,grid,figure,clf,stem,subplot)
from numpy import(abs,angle,pi,linspace)

deltat=.5
#n=6
#x=[0]*(2**n-ton)+[1]*2*ton+[0]*(2**n-ton)

t=linspace(-10,10,100)

def x(tlist):
    def xfun(t):
        if t>-deltat and t<deltat:
            return 1
        else:
            return 0
    return map(xfun,tlist)

figure(1)
clf()
subplot(2,1,1)
stem(t,x(t))
grid()

xt=fftshift(fft(x(t)))

subplot(4,1,3)
stem(abs(xt))
grid()
subplot(4,1,4)
stem(angle(xt)/pi)
grid()