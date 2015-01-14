from numpy.fft import (fft,fftshift)
from matplotlib.pyplot import (grid,figure,clf,stem,subplot,xlabel,xlim)
from numpy import(abs,angle,pi,linspace,cos,mod,concatenate)

def clean(a):
    '''Sets array values below tolerance to zero.'''
    tol=1e-13
    a.real[abs(a.real) < tol] = 0.0
    a.imag[abs(a.imag) < tol] = 0.0
    return a

def dftplot(t,x):
    '''Plots a sequence an its centered dft.'''
    
    fs=1/(t[1]-t[0])    
    fmin=-fs/2+mod(2*fs,2)/4.
    fmax=fs/2+mod(2*fs,2)/4.
    f=linspace(fmin,fmax,len(t))
    
    fig=figure(1)
    clf()
    subplot(2,1,1)
    stem(t,x)
    grid()
    xlabel('time (s)')
    
    x=x[:-1]
    xt=clean(fft(x))
    xt=fftshift(xt)
    xt=concatenate([xt,[xt[0]]])
    
    subplot(4,1,3)
    stem(f,abs(xt))
    #axvline(10)
    grid()
    xlim([fmin,fmax])
    subplot(4,1,4)
    stem(f,angle(xt)/pi)
    grid()
    xlim([fmin,fmax])    
    xlabel('frequency (Hz)')
    return fig
    
tstart=-1
tend=1
fs=15
bign=(tend-tstart)*fs+1
t=linspace(tstart,tend,bign)
x=cos(2*pi*1.5*t)#+pi/2)
dftplot(t,x)




