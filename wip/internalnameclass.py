"""
framework for the creation of an internal
namespace to customize code
Mar 2015
"""

from __future__ import print_function
import numpy as np
#import matplotlib.pyplot as plt
import sigsys2015rots as rs

class Mystuff():
    '''internal namespace'''
    pi=3.14
my = Mystuff()

#checking access to different variables
print(np.pi)
print(my.pi)

def foo(arg):
    return np.sin(arg)
    
rs.plot(foo)

t=np.linspace(-10,10,10)
rs.plot(t,foo)

rs.plot(t,t)

rs.plot(foo,color='r')

rs.plot(foo,figure=2)

rs.plot(foo,figure=3)