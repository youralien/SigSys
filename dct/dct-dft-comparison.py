'''comparison graphs between the dct and dft
jan 2015'''
from matplotlib.pyplot import (figure,subplots_adjust,clf,subplot,
                               stem,ylabel,xlabel,setp,plot,axhline,axvline,
                               grid,ylim,xlim,show,fill_between,yscale,legend)
from numpy import (array,set_printoptions,linspace,arange,angle,pi,
                   imag,real,concatenate)
#from cmath import * #why does this even exist?
from scipy.fftpack import dct
from myfunctions import (piecew,mydft,mydct4,myidft,myidct4)
# numpy output options
set_printoptions(precision=2,suppress=True)



'''Perform transforms and inverses below
x is an original function
DFTx is a DFT
iDFTx is an inverse DFT
DCTx is a DCT
iDCTx is an inverse DCT
'''

#choose the first experiment or the second
STARTWITHX=True
if STARTWITHX:
    '''First experiment:
    specify x, compute and compare the DFT and DCT, 
    then compute and compare x, iDFT and iDCT'''
    ### Test sequences
    sequences=[
    # [0] constant; try *1,2,3,4...
    [1,1]*1,
    # [1] periodic sequence at 1/2 Hz; try *1,2,3,4...
    [3,1]*1,
    # [2] as before but without DC component
    [-1,1]*4,
    # [3] 1/4 Hz
    [-1,1,1,-1]*2,
    # [4] pretty!
    [3,1,-3,7,3,1],
    # [5]
    [-3.,1,3,-1,-2,1,3,-3],
    # [6] energy conservation!
    [0]*10+[20]+[0]*9,
    ]
    sequence=array(sequences[0],dtype=float)    
    
    N=len(sequence)
    x=piecew(sequence)
    DFTx=mydft(x,N)
    iDFTx=myidft(DFTx,N)
    DCTx=mydct4(x,N)
    iDCTx=myidct4(DCTx,N)
else:
    '''second experiment:
    specify the DCT, compute x, compute the DFT 
    and compare them '''
    N=2# HACK WARNING: make sure N matches the length of the sequence below!!!
    iDCTx=myidct4(piecew([1,0]),N)#  start with a basis fn of DCT
    DCTx=mydct4(iDCTx,N)# take the DCT to check
    DFTx=mydft(iDCTx,N)# take the DFT
    iDFTx=myidft(DFTx,N)# take the iDFT back to check
    x=myidft(DFTx,N)# ditto



### setting the range
tstart=-2*N; tend=3*N # extended range to see periodicity
t=linspace(tstart,tend,500) # continuous range
nrange=arange(N) #discrete range
nextended=arange(tstart,tend) # extended discrete range to see periodicity
dctrange=arange(0+.25,N/2+.25,.5) # shifted range to compare the DCT and DFT



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
#close('all')
myfig=figure(1,figsize=(16,8), dpi=120)
clf()
subplots_adjust(wspace=.15,hspace=0.03)

#choosing colors and styles
originalfmt=[(0,.7,0),'s',14,'original x']
dftfmt=[(0,0,.6),'d',16,'DFT']
dctfmt=[(.6,0,0),'o',8,'DCT']
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
mystem(nrange,x(nrange),originalfmt)
#myplot(t,x(t),originalcolor) # continuous piecewise function

# DFT reconstruction
mystem(nextended,real(iDFTx(nextended)),dftfmt)


# DCT reconstruction
mystem(nextended,iDCTx(nextended),dctfmt)
axvline(x=-.5,color=dctfmt[0],linewidth=1,linestyle='--')

ylabel('x')
myscalemax=max(concatenate([iDCTx(t),abs(iDFTx(t))]))
myscalemin=min(concatenate([iDCTx(t),abs(iDFTx(t))]))
ylim(1.05*myscalemin,1.05*myscalemax)
setp(ax,xticklabels=[])
xlim(tstart,tend)
legend(framealpha=0,fontsize=10,loc=3,markerscale=.5)



### Time domain plot with interpolation
ax=subplot(223)
myaxes()

# original
mystem(nrange,x(nrange),originalfmt)
#myplot(t,x(t),originalcolor) # continuous piecewise function

# DFT reconstruction
mystem(nextended,real(iDFTx(nextended)),dftfmt)
myplot(t,real(iDFTx(t)),dftfmt)
fill_between(t,0,imag(iDFTx(t)),alpha=0.1,
             color=dftfmt[0])

# DCT reconstruction
mystem(nextended,iDCTx(nextended),dctfmt)
myplot(t,iDCTx(t),dctfmt)
axvline(x=-.5,color=dctfmt[0],linewidth=1,linestyle='--')

ylabel('x')
myscalemax=max(concatenate([iDCTx(t),abs(iDFTx(t))]))
myscalemin=min(concatenate([iDCTx(t),abs(iDFTx(t))]))
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
dctanddft=concatenate([abs(DCTx(t)),abs(DFTx(t))])
myscalemax=max(dctanddft)
myscalemin=max(10**-2,min(dctanddft))
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