from numpy import (arange,pi,exp,cos,floor,mod)

def piecew(seq):
    """
    Convert sequences to periodic piecewise functions.
    """
    def myfun(t):
        try:
            return [seq[int(floor(i))] for i in mod(t+.5,len(seq))]
            # t+.5 shifts to the center, can be changed
        except TypeError:
            return seq[int(floor(mod(t,len(seq))))]
    return myfun

def mydft(x,N):
    """
    Function implementation of DFT:
    DFT{x[n]} = 1/N \sum_n x[n]e^{-j 2pi f n/N}, n=0...N
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
