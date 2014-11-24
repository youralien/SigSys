from matplotlib.pyplot import (close,figure,subplots_adjust,get_cmap,subplot,
                               stem,ylabel,xlabel,setp,plot,axhline,axvline,
                               grid,ylim,xlim,show,fill_between,yscale,legend)
from numpy import (array,set_printoptions,linspace,arange,angle,pi,exp,cos,
                   imag,real,sqrt,concatenate,floor,mod)
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
sequence=array(sequences[0],dtype=float)    


#convert sequences to periodic piecewise functions
def piecew(seq):
    def out(t):
        try:
            return [seq[int(floor(i))] for i in mod(t+.5,len(seq))]
            # t+.5 shifts to the center, can be changed
        except TypeError:
            return seq[int(floor(mod(t,len(seq))))]
    return out



### setting the range
N=len(sequence); tstart=-1*N; tend=2*N # extended range to see periodicity
t=linspace(tstart,tend,500) # continuous range
nrange=arange(N) #discrete range
nextended=arange(tstart,tend) # extended discrete range to see periodicity
dctrange=arange(0+.25,N/2+.25,.5) # shifted range to compare the DCT and DFT



### Implementation of DFT, DCT type 4 and inverses
def mydft(x,N):
    """
    Function implementation of DFT:
    DFT{x[n]} = 1/N \sum_n x[n]e^{-j 2pi f n/N}, n=0...N
    Requires a global N.
    """
    def myfun(f):
        mysum=0
        for n in arange(N):
            mysum+=x(n)*exp(-1j*2.*pi*f*n/N)/N
        return mysum
    return myfun

### unused list implementation
#mydft=zeros(N,dtype=complex)
#for k in arange(N):
#    for n in arange(N):
#        mydft[k]+=x[n]*exp(-1j*2.*pi*k*n/N)/N

def mydct4(x,N):
    """
    Function implementation of DCT type 4:
    DCT{x[n]} = 1/N \sum_n x[n]cos(pi(f+0.5)(n+0.5)/N), n=0...N
    """
    def myfun(f):
        mysum=0
        for n in arange(N):
            mysum+=x(n)/N*cos(pi*(f+.5)*(n+.5)/N)
        return mysum
    return myfun
    
### unused list implementation
#mydct4=zeros(N)
#for k in arange(N):
#    for n in arange(N):
#        mydct4[k]+=x[n]/N*cos(pi*(k+.5)*(n+.5)/N)

def myidft(fx,N):
    """
    Function implementation of inverse DFT:
    IDFT{ DFT{x[n]} } = \sum_k DFT{x[n]}(k) e^{j 2pi k t/N}, k=0...N
    """
    def myfun(t):
        mysum=0
        for k in arange(N):
            mysum+=fx(k)*exp(1j*2.*pi*k*t/N)
        return mysum
    return myfun

def myidct4(fx,N):
    """
    Function implementation of inverse DCT type 4:
    IDCT{ DCT{x[n]} } = 2\sum_k DCT{x[n]}(k) cos(pi(k+0.5)(t+0.5)/N), k=0...N
    """
    def myfun(t):
        mysum=0
        for k in arange(N):
            mysum+=2*fx(k)*cos(pi*(k+.5)*(t+.5)/N)
        return mysum
    return myfun

#perform transforms and inverses    
x=piecew(sequence)  
DFTx=mydft(x,N)
DFTinvDFTx=myidft(DFTx,N)
#compute the norm of the imaginary part of the DFT reconstruction
imx=imag(myidft(nextended,N))
print '\n||Im(IDFT{ DFT{x} })||=',sqrt(imx.dot(imx))

DCTx=mydct4(x,N)
DCTinvDCTx=myidct4(DCTx,N)



### Analysis
# Console output of DCT types 1-4 and DFT
# Note: with scaling, I'm looking for factors that conserve 
# signal energy across domains.
print 'type 1:         ',dct(sequence,type=1)
print 'type 1 (scaled):',dct(sequence,type=1)*2/(1*N)
print 'type 2:         ',dct(sequence,type=2)
print 'type 2 (scaled):',dct(sequence,type=2)*2/(1*N)
print 'type 3:         ',dct(sequence,type=3)
print 'type 3 (scaled):',dct(sequence,type=3)*2/(1*N)
print 'type 4:         ',array(map(abs,DCTx(nrange)))
print 'dft:            ',array(map(abs,(DFTx(nrange))))



### Plotting
close('all')
myfig=figure(1,figsize=(16,8), dpi=120)
subplots_adjust(wspace=.15,hspace=0.03)

originalfmt=[(0,.7,0),'s',14,'original x']
dftfmt=[(0,0,.6),'d',16,'DFT']
dctfmt=[(.6,0,0),'o',8,'DCT']
#originalfmt=[(0,.7,0),'s',14,'original x']
#dftfmt=[(0,.4,.4),'d',16,'DFT']
#dctfmt=[(.1,0,.4),'o',8,'DCT']
#spectral colors go from purple to red, 0 is black
#dftcolor=get_cmap('nipy_spectral')(.2)[:3]
psdcolor=(.5,.2,.2)
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

def mystem(x,y,myformat):
    markerline, stemlines, baseline = stem(x,y, '-.',label=myformat[3],
                                           markerfmt=myformat[1])
    setp(markerline,'markerfacecolor',myformat[0]+(.4,),'markersize',myformat[2],
         'markeredgecolor',myformat[0]+(1,),'markeredgewidth',1)
    setp(stemlines,'color',myformat[0]+(.5,))
    setp(baseline,'color',myformat[0]+(0,))
    

def myplot(x,y,myformat):
    lines=plot(x,y,'-')
    setp(lines,'color',myformat[0],'linewidth',2,'alpha',.7)

### Time domain plot
ax=subplot(221)
myaxes()

# original
mystem(nrange,sequence,originalfmt)
#myplot(t,x(t),originalcolor) # plotting continuous piecewise function

# DFT reconstruction
mystem(nextended,real(DFTinvDFTx(nextended)),dftfmt)


# DCT reconstruction
mystem(nextended,DCTinvDCTx(nextended),dctfmt)
axvline(x=-.5,color=dctfmt[0],linewidth=1,linestyle='--')

ylabel('x')
myscalemax=max(concatenate([DCTinvDCTx(t),abs(DFTinvDFTx(t))]))
myscalemin=min(concatenate([DCTinvDCTx(t),abs(DFTinvDFTx(t))]))
ylim(1.05*myscalemin,1.05*myscalemax)
setp(ax,xticklabels=[])
xlim(tstart,tend)
legend(framealpha=0,fontsize=10,loc=3,markerscale=.5)



### Time domain plot with interpolation
ax=subplot(223)
myaxes()

# original
mystem(nrange,sequence,originalfmt)
#myplot(t,x(t),originalcolor) # continuous piecewise function

# DFT reconstruction
mystem(nextended,real(DFTinvDFTx(nextended)),dftfmt)
myplot(t,real(DFTinvDFTx(t)),dftfmt)
#myplot(t,imag(myidft(t)),[min(1,c+.5) for c in dftcolor]) #imaginary line
fill_between(t,real(DFTinvDFTx(t)),imag(DFTinvDFTx(t)),alpha=0.1,
             color=dftfmt[0])

# DCT reconstruction
mystem(nextended,DCTinvDCTx(nextended),dctfmt)
myplot(t,DCTinvDCTx(t),dctfmt)
axvline(x=-.5,color=dctfmt[0],linewidth=1,linestyle='--')

ylabel('x')
myscalemax=max(concatenate([DCTinvDCTx(t),abs(DFTinvDFTx(t))]))
myscalemin=min(concatenate([DCTinvDCTx(t),abs(DFTinvDFTx(t))]))
ylim(1.05*myscalemin,1.05*myscalemax)
xlabel('time (s)')
xlim(tstart,tend)
legend(framealpha=0,fontsize=10,loc=3,markerscale=.5)



### Magnitude
ax=subplot(222)  
myaxes()

mystem(nrange,abs(DFTx(nrange)),dftfmt)
myplot(t,abs(DFTx(t)),dftfmt)
mystem(dctrange,abs(DCTx(nrange)),dctfmt)
myplot(t/2+.25,abs(DCTx(t)),dctfmt)

ylabel('magnitude of FT{x}')
yscale('log')
myscalemax=max(concatenate([abs(DCTx(t)),abs(DFTx(t))]))
myscalemin=min(abs(DCTx(nextended)))
ylim(.8*myscalemin,1.2*myscalemax)
xlim(-1,N)
legend(framealpha=0,fontsize=10,loc=3,markerscale=.5)
setp(ax,xticklabels=[])
#freqaxis() # use the axes from the subplot below

#secondary scale for PSD
ax2=ax.twinx()
ylabel('PSD of x',color=psdcolor)
for label in ax2.get_yticklabels():
    label.set_color(psdcolor)
yscale('log')
ylim((.8*myscalemin)**2,(1.2*myscalemax)**2)
grid(which='both',color=psdcolor,linestyle='-',linewidth=.2)




### Angle
ax=subplot(224)
myaxes()

#mystem(nrange,map(lambda aux:angle(mydft(x)(aux))/pi,nrange),dftcolor)
mystem(nrange,angle(DFTx(nrange))/pi,dftfmt)
myplot(t,angle(DFTx(t))/pi,dftfmt)
mystem(dctrange,angle(DCTx(nrange))/pi,dctfmt)
myplot(t/2+.25,angle(DCTx(t))/pi,dctfmt)

ylabel('angle/$\pi$ of FT{x}')
ylim(-1.2,1.2)
freqaxis()
legend(framealpha=0,fontsize=10,loc=3,markerscale=.5)



### unused PSD plot (not necessary since magnitude can simply be stretched)
#ax=subplot(223)
#myaxes()
#  
#def dftpsd(x):return abs(mydft(x))**2
#mystem(nrange,map(dftpsd,nrange),dftcolor)
#myplot(t,map(dftpsd,t),dftcolor)
#
#def dctpsd(x):return abs(mydct4(x))**2
#mystem(dctrange,map(dctpsd,nrange),dctcolor)
#myplot(t/2+.25,map(dctpsd,t),dctcolor)
#
#ylabel('PSD of x')
#yscale('log')
#myscalemax=max(concatenate([abs(dctpsd(t)),abs(dftpsd(t))]))
#myscalemin=min(abs(dctpsd(nextended)))
#ylim(.8*myscalemin,1.2*myscalemax)
#freqaxis()