from matplotlib.pyplot import *
from numpy import *
from cmath import *
from scipy.fftpack import dct

set_printoptions(precision=2,suppress=True)


#periodic sequence at 1/2 Hz
x=array([3.,1.])
x=array([3.,1.,3.,1.])
x=array([3.,1.,3.,1.,3.,1.])
x=array([3.,1.,3.,1.,3.,1.,3.,1.])
x=array([-1.,1.,-1.,1.,-1.,1.,-1.,1.])
x=array([-1.,1.,1.,-1.,-1.,1.,1.,-1.])


x=array([3.,1.,-3.,7.,3.,1.])#pretty!

x=array([-3.,1,3,-1,-2,1,3,-3])

#x=array([0.,0,0,0,0,0,0,0,0,0,20,0,0,0,0,0,0,0,0,0])#energy conservation!



N=len(x)

tstart=-1*N
tend=2*N
t=linspace(tstart,tend,500)
nrange=arange(tstart,tend)

#ftx=zeros(N,dtype=complex)
#for k in arange(N):
#    for n in arange(N):
#        ftx[k]+=x[n]*exp(-1j*2.*pi*k*n/N)/N

def ftx(f):
    ft=0
    for n in arange(N):
        ft+=x[n]*exp(-1j*2.*pi*f*n/N)/N
    return ft

#dct4x=zeros(N)
#for k in arange(N):
#    for n in arange(N):
#        dct4x[k]+=x[n]/N*cos(pi*(k+.5)*(n+.5)/N)

def dct4x(f):
    ft=0
    for n in arange(N):
        ft+=x[n]/N*cos(pi*(f+.5)*(n+.5)/N)
    return ft
    
dctrange=arange(0+.25,N/2+.25,.5)

print dct(x,type=1),'type 1'
print dct(x,type=1)*2/(1*N),'type 1 (scaled)'
print dct(x,type=2),'type 2:'
print dct(x,type=2)*2/(1*N),'type 2 (scaled)'
print dct(x,type=3),'type 3:'
print dct(x,type=3)*2/(1*N),'type 3 (scaled)'
print array(map(abs,map(dct4x,arange(N)))),'type 4'
print array(map(abs,map(ftx,arange(N)))),'dft'

def xft(t):
    x=0
    for k in arange(N):
        x+=ftx(k)*exp(1j*2.*pi*k*t/N)
    return x

def xdct4(t):
    x=0
    for k in arange(N):
        x+=2*dct4x(k)*cos(pi*(k+.5)*(t+.5)/N)
    return x




################# Plotting
close('all')
myfig=figure(1,figsize=(16,8), dpi=120)
subplots_adjust(wspace=.15,hspace=.3)
mycolor=get_cmap('nipy_spectral')(.1)[:3]





ax=subplot(321)  

markerline, stemlines, baseline = stem(x, '-.')
ylabel('x (original)')
xlabel('time (s)')

setp(markerline, 'markerfacecolor', 'b')
setp(stemlines,'color','b')
axhline(color='k')
axvline(color='k')
grid(which='both')
ylim(min(-1,min(x)-1),max(x)+1)
xlim(tstart,tend)








ax=subplot(322)  

markerline, stemlines, baseline = stem(map(abs,map(ftx,arange(N))), '-.')
plot(t,map(abs,map(ftx,t)),'b')
markerline2, stemlines2, baseline2 = stem(dctrange,map(abs,map(dct4x,arange(0,N))), '-.')
plot(t/2+.25,map(abs,map(dct4x,t)),'r')
ylabel('magnitude of FT{x}')
xlabel('frequency (Hz)')
ax.set_xticks(arange(-1,N+1))
ax.set_xticklabels([r'$\frac{'+str(n)+'}{'+str(N)+'}$' for n in arange(-1,N+1)])

setp(markerline, 'markerfacecolor', 'b')
setp(stemlines,'color','b')
setp(markerline2, 'markerfacecolor', 'r')
setp(stemlines2,'color','r')
axhline(color='k')
axvline(color='k')
grid(which='both')
ylim(-.2,max(map(dctpsd,t)+map(ftpsd,t))+.2)
xlim(-1,N)


ax=subplot(324)

markerline, stemlines, baseline = stem(map(lambda x:angle(ftx(x))/pi,arange(N)), '-.')
plot(t,map(lambda x:angle(ftx(x))/pi,t),'b')
markerline2, stemlines2, baseline2 = stem(dctrange,map(lambda x:angle(dct4x(x))/pi,arange(0,N)), '-.')
plot(t/2+.25,map(lambda x:angle(dct4x(x))/pi,t),'r')
ylabel('angle/$\pi$ of FT{x}')
xlabel('frequency (Hz)')
ax.set_xticks(arange(-1,N+1))
ax.set_xticklabels([r'$\frac{'+str(n)+'}{'+str(N)+'}$' for n in arange(-1,N+1)])

setp(markerline, 'markerfacecolor', 'b')
setp(stemlines,'color','b')
setp(markerline2, 'markerfacecolor', 'r')
setp(stemlines2,'color','r')
axhline(color='k')
axvline(color='k')
grid(which='both')
ylim(-1.2,1.2)
xlim(-1,N)


ax=subplot(313)  
ftpsd=lambda x:abs(ftx(x))**2
markerline, stemlines, baseline = stem(map(ftpsd,arange(N)), '-.')
plot(t,map(ftpsd,t),'b')
dctpsd=lambda x:abs(dct4x(x))**2
markerline2, stemlines2, baseline2 = stem(dctrange,map(dctpsd,arange(0,N)), '-.')
plot(t/2+.25,map(dctpsd,t),'r')
ylabel('PSD of x')
xlabel('frequency (Hz)')
ax.set_xticks(arange(-1,N+1))
ax.set_xticklabels([r'$\frac{'+str(n)+'}{'+str(N)+'}$' for n in arange(-1,N+1)])

setp(markerline, 'markerfacecolor', 'b')
setp(stemlines,'color','b')
setp(markerline2, 'markerfacecolor', 'r')
setp(stemlines2,'color','r')
axhline(color='k')
axvline(color='k')
grid(which='both')
ylim(-.2,max(map(dctpsd,t)+map(ftpsd,t))+.2)
xlim(-1,N)



ax=subplot(323)
markerline, stemlines, baseline = stem(nrange,map(xft,nrange), '-.')
plot(t,map(xft,t),'b')
markerline2, stemlines2, baseline2 = stem(nrange,map(xdct4,nrange), '-.')
plot(t,map(xdct4,t),'r')
ylabel('x (reconstruted)')
xlabel('time (s)')

setp(markerline, 'markerfacecolor', 'b')
setp(stemlines,'color','b')
setp(markerline2, 'markerfacecolor', 'r')
setp(stemlines2,'color','r')
axhline(color='k')
axvline(color='k')
axvline(x=-.5,color='r',linewidth=1,linestyle='--')
grid(which='both')
ylim(min(-1,min(map(xdct4,t))-1),max(map(xdct4,t))+1)
xlim(tstart,tend)


show()

