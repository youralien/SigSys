from numpy.fft import (fft,fftshift)
from matplotlib.pyplot import (grid,figure,clf,stem,subplot,xlabel,xlim,plot,
                               setp,ylabel,ylim)
from numpy import(abs,angle,pi,linspace,cos,mod,concatenate,array)

def clean(a):
    '''Sets array values below tolerance to zero.'''
    tol=1e-13
    a.real[abs(a.real) < tol] = 0.0
    a.imag[abs(a.imag) < tol] = 0.0
    return a

def dftplot(t,x):
    '''Plots a sequence an its centered dft.'''
    
    fs=1/(t[1]-t[0])    
    fmin=-fs/2+mod(fs,2)/2.
    fmax=fs/2+mod(fs,2)/2.
    f=linspace(fmin,fmax,len(t))
    
    def mystem(t,x):
        markerline, stemlines, baseline = stem(t,x)#, markerfmt="")
        setp(baseline, 'color','k', 'linewidth', 1)        
        #setp(markerline, 'color',color, 'linewidth', 1)       
        
    fig=figure(1)
    clf()
    subplot(2,1,1)
    mystem(t[:-1],x[:-1])
    mystem([t[-1]],[x[-1]])
    #plot(t,x,':',lw=.5)
    grid()
    ylabel('signal')
    xlabel('time (s)')
    
    x=x[:-1]
    xt=clean(fft(x)/(len(t)-1))
    xt=fftshift(xt)
    xt=concatenate([xt,[xt[0]]])
    
    subplot(4,1,3)
    mystem(f,abs(xt))
    #plot(f,abs(xt),':',lw=.5)
    #axvline(10)
    grid()
    ylabel('magnitude')
    xlim([1.1*fmin,1.1*fmax])
    ylim([0,1.1*max(abs(xt))])
    subplot(4,1,4)
    mystem(f,angle(xt)/pi)
    #plot(f,angle(xt)/pi,':',lw=.5)
    grid()
    ylabel('phase/pi')
    xlim([1.1*fmin,1.1*fmax])
    ylim([-1,1])    
    xlabel('frequency (Hz)')
    return fig
    
tstart=0
tend=1
fs=2
bign=(tend-tstart)*fs+1
t=linspace(tstart,tend,bign)
x=cos(2*pi*1*t+0)#+pi/2)
dftplot(t,x)




