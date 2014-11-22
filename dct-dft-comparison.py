from matplotlib.pyplot import (close,figure,subplots_adjust,get_cmap,subplot,
                               stem,ylabel,xlabel,setp,plot,axhline,axvline,
                               grid,ylim,xlim,show,fill_between,yscale)
from numpy import (array,set_printoptions,linspace,arange,angle,pi,exp,cos,
                   imag,real,sqrt,concatenate)
#from cmath import * #why does this even exist?
from scipy.fftpack import dct

# numpy output options
set_printoptions(precision=2,suppress=True)



### Test sequences
sequences=[
# [0] periodic sequence at 1/2 Hz; try *1,2,3,4...
[3,1]*1,
# [1] without DC component
[-1,1]*4,
# [2] 1/4 Hz
[-1,1,1,-1]*2,
# [3] pretty!
[3,1,-3,7,3,1],
# [4]
[-3.,1,3,-1,-2,1,3,-3],
# [5] energy conservation!
[0]*10+[20]+[0]*9,
]
x=array(sequences[0],dtype=float)



### setting the range
N=len(x); tstart=-1*N; tend=2*N
t=linspace(tstart,tend,500) # continuous range
nrange=arange(tstart,tend) # discrete range
dctrange=arange(0+.25,N/2+.25,.5) # shifted range to compare the DCT and DFT



### Implementation of DFT, DCT type 4 and inverses
def mydft(f):
    """
    Function implementation of DFT:
    DFT{x[n]} = 1/N \sum_n x[n]e^{-j 2pi f n/N}, n=0...N
    Requires a global x[n].
    """
    mysum=0
    for n in arange(N):
        mysum+=x[n]*exp(-1j*2.*pi*f*n/N)/N
    return mysum
    
### unused list implementation
#mydft=zeros(N,dtype=complex)
#for k in arange(N):
#    for n in arange(N):
#        mydft[k]+=x[n]*exp(-1j*2.*pi*k*n/N)/N

def mydct4(f):
    """
    Function implementation of DCT type 4:
    DCT{x[n]} = 1/N \sum_n x[n]cos(pi(f+0.5)(n+0.5)/N), n=0...N
    Requires a global x[n].
    """
    mysum=0
    for n in arange(N):
        mysum+=x[n]/N*cos(pi*(f+.5)*(n+.5)/N)
    return mysum
    
### unused list implementation
#mydct4=zeros(N)
#for k in arange(N):
#    for n in arange(N):
#        mydct4[k]+=x[n]/N*cos(pi*(k+.5)*(n+.5)/N)

def myidft(t):
    """
    Function implementation of inverse DFT:
    IDFT{ DFT{x[n]} } = \sum_k DFT{x[n]}(k) e^{j 2pi k t/N}, k=0...N
    Reqquires a global mydft(f) function.
    """
    mysum=0
    for k in arange(N):
        mysum+=mydft(k)*exp(1j*2.*pi*k*t/N)
    return mysum

def myidct4(t):
    """
    Function implementation of inverse DCT type 4:
    IDCT{ DCT{x[n]} } = 2\sum_k DCT{x[n]}(k) cos(pi(k+0.5)(t+0.5)/N), k=0...N
    Requires a global mydct(f) function.
    """
    mysum=0
    for k in arange(N):
        mysum+=2*mydct4(k)*cos(pi*(k+.5)*(t+.5)/N)
    return mysum



### Analysis
# Console output of DCT types 1-4 and DFT
# Note: with scaling, I'm looking for factors that conserve 
# signal energy across domains.
print 'type 1:         ',dct(x,type=1)
print 'type 1 (scaled):',dct(x,type=1)*2/(1*N)
print 'type 2:         ',dct(x,type=2)
print 'type 2 (scaled):',dct(x,type=2)*2/(1*N)
print 'type 3:         ',dct(x,type=3)
print 'type 3 (scaled):',dct(x,type=3)*2/(1*N)
print 'type 4:         ',array(map(abs,map(mydct4,arange(N))))
print 'dft:            ',array(map(abs,map(mydft,arange(N))))



### Plotting
close('all')
myfig=figure(1,figsize=(16,8), dpi=120)
subplots_adjust(wspace=.15,hspace=0.02)
#spectral colors go from purple to red, 0 is black
dftcolor=get_cmap('nipy_spectral')(.2)[:3]
dctcolor=get_cmap('nipy_spectral')(.9)[:3]
originalcolor=get_cmap('nipy_spectral')(.6)[:3]
psdcolor=get_cmap('nipy_spectral')(.1)[:3]
show()

def myaxes():
    axhline(color='k')
    axvline(color='k')
    grid(which='both')
    
def freqaxis():
    xlabel('frequency (Hz)')
    xlim(-1,N)
    ax.set_xticks(arange(-1,N+1))
    ax.set_xticklabels([r'$\frac{'+str(n)+'}{'+str(N)+'}$' 
                        for n in arange(-1,N+1)])

def mystem(x,y,mycolor):
    markerline, stemlines, baseline = stem(x,y, '-.')
    setp(markerline,'markerfacecolor',mycolor+(.6,),'markersize',12)
    setp(stemlines,'color',mycolor+(.5,))
    setp(baseline,'color',mycolor+(.5,))

def myplot(x,y,mycolor):
    lines=plot(x,y,'-')
    setp(lines,'color',mycolor)



### Time domain plot
ax=subplot(121)
myaxes()

# DFT reconstruction
imx=imag(myidft(nrange))
print '\n||Im(IDFT{ DFT{x} })||=',sqrt(imx.dot(imx))
mystem(nrange,real(myidft(nrange)),dftcolor)
myplot(t,real(myidft(t)),dftcolor)
#myplot(t,imag(myidft(t)),[min(1,c+.5) for c in dftcolor])
fill_between(t,real(myidft(t)),imag(myidft(t)),alpha=0.1,color=dftcolor)

# DCT reconstruction
mystem(nrange,myidct4(nrange), dctcolor)
myplot(t,myidct4(t),dctcolor)
axvline(x=-.5,color=dctcolor,linewidth=1,linestyle='--')

# original
mystem(arange(N),x,originalcolor)

ylabel('x')
myscalemax=max(concatenate([myidct4(t),abs(myidft(t))]))
myscalemin=min(concatenate([myidct4(t),abs(myidft(t))]))
ylim(1.05*myscalemin,1.05*myscalemax)
xlabel('time (s)')
xlim(tstart,tend)



### Magnitude
ax=subplot(222)  
myaxes()

mystem(arange(N),map(abs,map(mydft,arange(N))),dftcolor)
myplot(t,map(abs,map(mydft,t)),dftcolor)
mystem(dctrange,map(abs,map(mydct4,arange(0,N))),dctcolor)
myplot(t/2+.25,map(abs,map(mydct4,t)),dctcolor)

ylabel('magnitude of FT{x}')
yscale('log')
myscalemax=max(concatenate([abs(mydct4(t)),abs(mydft(t))]))
myscalemin=min(abs(mydct4(nrange)))
ylim(.8*myscalemin,1.2*myscalemax)
xlim(-1,N)
setp(ax,xticklabels=[])
#freqaxis() # use the axes from the subplot below

#secondary scale for PSD
ax2=ax.twinx()
ylabel('PSD of x',color=psdcolor)
for label in ax2.get_yticklabels():
    label.set_color(psdcolor)
yscale('log')
ylim((.8*myscalemin)**2,(1.2*myscalemax)**2)
grid(which='both',color=psdcolor,linestyle='-',linewidth=.1)



### Angle
ax=subplot(224)
myaxes()

mystem(arange(N),map(lambda x:angle(mydft(x))/pi,arange(N)),dftcolor)
myplot(t,map(lambda x:angle(mydft(x))/pi,t),dftcolor)
mystem(dctrange,map(lambda x:angle(mydct4(x))/pi,arange(0,N)),dctcolor)
myplot(t/2+.25,map(lambda x:angle(mydct4(x))/pi,t),dctcolor)

ylabel('angle/$\pi$ of FT{x}')
ylim(-1.2,1.2)
freqaxis()



### unused PSD plot (not necessary since magnitude can simply be stretched)
#ax=subplot(223)
#myaxes()
#  
#def dftpsd(x):return abs(mydft(x))**2
#mystem(arange(N),map(dftpsd,arange(N)),dftcolor)
#myplot(t,map(dftpsd,t),dftcolor)
#
#def dctpsd(x):return abs(mydct4(x))**2
#mystem(dctrange,map(dctpsd,arange(0,N)),dctcolor)
#myplot(t/2+.25,map(dctpsd,t),dctcolor)
#
#ylabel('PSD of x')
#yscale('log')
#myscalemax=max(concatenate([abs(dctpsd(t)),abs(dftpsd(t))]))
#myscalemin=min(abs(dctpsd(nrange)))
#ylim(.8*myscalemin,1.2*myscalemax)
#freqaxis()