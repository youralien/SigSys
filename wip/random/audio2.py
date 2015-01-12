import numpy as np
from matplotlib.pyplot import plot
from scipy.io import wavfile

def show_info(a):
    print "shape:", a.shape
    print "dtype:", a.dtype
    print "min, max:", a.min(), a.max()
    print

rate, data = wavfile.read('computermercy2.wav')

show_info(data)
plot(data)
## Take the sine of each element in `data`.
## The np.sin function is "vectorized", so there is no need
## for a Python loop here.
#sindata = np.sin(data)
#
#show_info("sindata", sindata)
#
## Scale up the values to 16 bit integer range and round
## the value.
#scaled = np.round(32767*sindata)
#
#show_info("scaled", scaled)
#
## Cast `scaled` to an array with a 16 bit signed integer data type.
#newdata = data.astype(np.int16)
#
#show_info("newdata", newdata)
#
## Write the data to 'newname.wav'
wavfile.write('newname.wav', rate, data)
