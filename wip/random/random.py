from numpy import (linspace,sqrt)
from numpy.random import (normal,uniform)
from matplotlib.pyplot import (plot,xlabel,ylabel,grid,show,
                               clf,figure,ylim)
import winsound

#winsound.PlaySound('computer_mercy.wav',winsound.SND_FILENAME)

data=open('c:/windows/media/Chimes.wav',"rb").read()
print data
#winsound.PlaySound(data, winsound.SND_MEMORY)

N=1e3
t=linspace(0,10,N)
xnor = normal(size=N)
xuni = uniform(-sqrt(3),sqrt(3),size=N)


### Plotting
myfig=figure(1,figsize=(16,8), dpi=120)
clf()
#subplots_adjust(wspace=.15,hspace=0.03)

plot(t,xnor)
plot(t,xuni)
xlabel('time')
ylabel('x')
ylim([-4,4])
grid()
show()


