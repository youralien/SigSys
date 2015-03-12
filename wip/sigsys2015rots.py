"""
Olin College code for 
sigsys 2015 episode 3: Revenge of the Sith
Mar 2015
jomm@olin.edu
"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

class MyError(TypeError):
    pass

def plot(*args,**kwargs):
    '''internal plot function
    first or second argument must be a function'''
    #print('\nplot is called')
    plt.close(1)
    if 'figure' in kwargs:
        plt.figure(kwargs['figure'])
    else:
        plt.figure(1)
    plt.grid(True)
    if callable(args[0]):
        #print('first arg is callable')
        #not smart: use an arbitrary range
        t=np.linspace(1,10,100)
        #only passing kwargs
        plt.plot(t,args[0](t),**kwargs)
    else:
        try:
            #print('first arg is not callable')
            #only passing kwargs
            plt.plot(args[0],args[1](args[0]),**kwargs)      
        except TypeError:
            raise MyException('bleh')
            print('\nfirst or second argument must be a function')
            print('args[0]:',type(args[0]))
            print('args[1]:',type(args[1]))
            plt.close(1)
            #raise TypeError
