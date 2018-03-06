# Using Delauney triangulation to look at how number of particles dropped into cylinder affects interparticle separation
import numpy
import pickle
import math
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import itertools
import constants as const
import matplotlib



filehandler = open("/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/crystalpositions2,5K.obj",
                   'rb')  ##2k particles 5.5hrs to run 2500 iterations just for settling down
xinit = pickle.load(filehandler)
yinit = pickle.load(filehandler)
zinit = pickle.load(filehandler)
filehandler.close()
initpositions = numpy.array([[i,j] for i,j in zip(xinit, yinit)])
triangulate = Delaunay(initpositions)

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 15}
matplotlib.rc('font', **font)
# Plot the triangulation
plt.triplot(initpositions[:,0], initpositions[:,1], triangulate.simplices.copy(),color='g')
plt.scatter(initpositions[:,0], initpositions[:,1], color='g',s=15)
plt.title("Triangulation of a section of a \n crystal composed of 2.5K dust particles")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.axis("equal")
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlim([-0.0015,0.0015])
plt.ylim([-0.0015,0.0015])
from pylab import rcParams
rcParams['figure.figsize'] = 9, 5
from pylab import rcParams
plt.savefig('triangulate', format='eps')

plt.show()

#Plot distribution of interparticle separation distances


edgelengths = []
pairlist = []
pairlistfinal = []

for i in numpy.arange(len(triangulate.simplices)):
    pairlist.extend(itertools.combinations(triangulate.simplices[i], 2))

pairlistfinal.extend(map(tuple, set(map(frozenset, pairlist))))

for j in numpy.arange(len(pairlistfinal)):
    indices = pairlistfinal[j]
    ind1 = initpositions[indices[0]]
    ind2 = initpositions[indices[1]]
    edgelengths.append(math.hypot(ind1[0] - ind2[0], ind1[1] - ind2[1]))

binwidth = const.lambdaD/5
plt.hist(edgelengths,bins=numpy.arange(min(edgelengths), max(edgelengths) + binwidth, binwidth))
plt.ylabel("Frequency")
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.xlabel("inter-dust separation distance (m)")
plt.title("Histogram of inter-dust separation distances of crystal")
plt.xlim([0.1*10**(-3),0.6*10**(-3)])
from pylab import rcParams
rcParams['figure.figsize'] = 9, 5
plt.savefig('histogram', format='eps')
plt.show()