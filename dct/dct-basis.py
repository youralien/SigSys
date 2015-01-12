'''plot of DCT basis functions
jan 2015'''
from matplotlib.pyplot import (clf,figure,subplots_adjust,get_cmap,subplot,
                               stem,ylabel,xlabel,setp,plot,axhline,axvline,
                               grid,ylim,xlim,show)
from numpy import (set_printoptions,linspace,arange)
from myfunctions import (piecew,mydft,myidft,myidct4)
# numpy output options
set_printoptions(precision=2,suppress=True)



### setting the range
N=4#N=len(sequence)
tstart=-1*N; tend=2*N # extended range to see periodicity
t=linspace(tstart,tend,500) # continuous range
nrange=arange(N) #discrete range
nextended=arange(tstart,tend) # extended discrete range to see periodicity
dctrange=arange(0+.25,N/2+.25,.5) # shifted range to compare the DCT and DFT

###perform transforms and inverses
#x is an original function
#DFTx is a DFT
#iDFTx is an inverse DFT
#DCTx is a DCT
#iDCTx is an inverse DCT

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



### Plotting
myfig=figure(1,figsize=(16,8), dpi=120)
clf()
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