# Using Delauney triangulation to look at how number of particles dropped into cylinder affects interparticle separation
import numpy
import pickle
import math
import matplotlib.pyplot as plt

from scipy.spatial import Delaunay


filehandler = open("/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/crystalpositions2,5K.obj",
                   'rb')  ##2k particles 5.5hrs to run 2500 iterations just for settling down
xinit = pickle.load(filehandler)
yinit = pickle.load(filehandler)
zinit = pickle.load(filehandler)
filehandler.close()
initpositions = numpy.array([[i,j] for i,j in zip(xinit, yinit)])

triangulate = Delaunay(initpositions)

#Plot distribution of interparticle separation distances

plt.triplot(initpositions[:,0], initpositions[:,1], triangulate.simplices.copy())
plt.plot(initpositions[:,0], initpositions[:,1], 'o')
plt.show()