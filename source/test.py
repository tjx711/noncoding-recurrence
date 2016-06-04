#!python

import numpy
from matplotlib import pyplot

x = numpy.arange(10)
y = numpy.array([5,3,4,2,7,5,4,6,3,2])

fig = pyplot.figure()
ax = fig.add_subplot(111)
ax.set_ylim(0,10)

pyplot.plot(y)
for i,j in zip(x,y):
    ax.annotate(str(j),xy=(i,j))


pyplot.show()
