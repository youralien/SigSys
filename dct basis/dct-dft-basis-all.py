from matplotlib.pyplot import (close,figure,subplots_adjust,get_cmap,subplot,
                               stem,ylabel,xlabel,setp,plot,axhline,axvline,
                               grid,ylim,xlim,show,fill_between,yscale,legend)
from numpy import (array,set_printoptions,linspace,arange,angle,pi,exp,cos,
                   imag,real,sqrt,concatenate,floor,mod)
# numpy output options
set_printoptions(precision=2,suppress=True)



#convert sequences to periodic piecewise functions
def piecew(seq):
    def out(t):
        try:
            return [seq[int(floor(i))] for i in mod(t+.5,len(seq))]
            # t+.5 shifts to the center, can be changed
        except TypeError:
            return seq[int(floor(mod(t,len(seq))))]
    return out



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



### setting the range
N=4#N=len(sequence)
tstart=-1*N; tend=2*N # extended range to see periodicity
t=linspace(tstart,tend,500) # continuous range
nrange=arange(N) #discrete range
nextended=arange(tstart,tend) # extended discrete range to see periodicity
dctrange=arange(0+.25,N/2+.25,.5) # shifted range to compare the DCT and DFT



iDCTx0=myidct4(piecew([1,0,0,0]),N)#  basis fn of DCT
iDCTx1=myidct4(piecew([0,1,0,0]),N)#  basis fn of DCT
iDCTx2=myidct4(piecew([0,0,1,0]),N)#  basis fn of DCT
iDCTx3=myidct4(piecew([0,0,0,1]),N)#  basis fn of DCT
DFTx0=mydft(iDCTx0,N)# take the DFT
DFTx1=mydft(iDCTx1,N)# take the DFT
DFTx2=mydft(iDCTx2,N)# take the DFT
DFTx3=mydft(iDCTx3,N)# take the DFT
x0=myidft(DFTx0,N)# check
x1=myidft(DFTx1,N)# check
x2=myidft(DFTx2,N)# check
x3=myidft(DFTx3,N)# check



###perform transforms and inverses
#x is an original function
#DFTx is a DFT
#iDFTx is an inverse DFT
#DCTx is a DCT
#iDCTx is an inverse DCT
#canonical order
#x=piecew([2,1.6,1.1,.4])
#DFTx=mydft(x,N)
#iDFTx=myidft(DFTx,N)
#DCTx=mydct4(x,N)
#iDCTx=myidct4(DCTx,N)



### Plotting
close('all')
myfig=figure(1,figsize=(16,8), dpi=120)
subplots_adjust(wspace=.15,hspace=0.03)

originalfmt=[(0,.5,0),'s',16,'original x']
mycmap='nipy_spectral_r';start=.1
dctfmt0=[(.8,0,0),'o',10,'DCT']
dctfmt1=[get_cmap(mycmap)(start+.1)[:3],'o',10,'DCT']
dctfmt2=[get_cmap(mycmap)(start+.7)[:3],'o',10,'DCT']
dctfmt3=[get_cmap(mycmap)(start+.8)[:3],'o',10,'DCT']
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
    setp(markerline,'markerfacecolor',myformat[0]+(.3,),'markersize',myformat[2],
         'markeredgecolor',myformat[0]+(1,),'markeredgewidth',1)
    setp(stemlines,'color',myformat[0]+(.5,))
    setp(baseline,'color',myformat[0]+(0,))
    

def myplot(x,y,myformat):
    lines=plot(x,y,'-')
    setp(lines,'color',myformat[0],'linewidth',4,'alpha',.99)



### Time domain plot with interpolation
ax=subplot(111)
myaxes()

# original
mystem(nrange,x0(nrange),originalfmt)
mystem(nrange,x1(nrange),originalfmt)
mystem(nrange,x2(nrange),originalfmt)
mystem(nrange,x3(nrange),originalfmt)

# DCT reconstruction
mystem(nextended,iDCTx0(nextended),dctfmt0)
mystem(nextended,iDCTx1(nextended),dctfmt1)
mystem(nextended,iDCTx2(nextended),dctfmt2)
mystem(nextended,iDCTx3(nextended),dctfmt3)
myplot(t,iDCTx0(t),dctfmt0)
myplot(t,iDCTx1(t),dctfmt1)
myplot(t,iDCTx2(t),dctfmt2)
myplot(t,iDCTx3(t),dctfmt3)
axvline(x=-.5,color=dctfmt0[0],linewidth=1,linestyle='--')

ylabel('x')
myscalemax=max(iDCTx0(t))
myscalemin=min(iDCTx0(t))
ylim(1.05*myscalemin,1.05*myscalemax)
xlabel('time (s)')
xlim(tstart,tend)

