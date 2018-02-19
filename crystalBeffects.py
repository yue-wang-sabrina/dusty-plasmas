# This file is mainly for looking at B field effects of crystals that are already formed.

import numpy
import matplotlib.pyplot as plt
import math
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
import pandas as pd
import mpl_toolkits.mplot3d.axes3d as p3
import itertools
import time
from scipy import special
import pickle
from tqdm import tqdm
from numpy import isclose
from scipy import interpolate


##Generate particles in their equilibrium position
filehandler = open("crystalpositions2K.obj", 'rb')  ##2k particles 5.5hrs to run 2500 iterations just for settling down
xinit = pickle.load(filehandler)
yinit = pickle.load(filehandler)
zinit = pickle.load(filehandler)
filehandler.close()
initpositions = [[i, j, k] for i, j, k in zip(xinit, yinit, zinit)]

## Quick plot to view initial configuration
# plt.figure()
# plt.scatter(xinit,yinit)
# plt.xlim((min(xinit),max(xinit)))
# plt.ylim((min(yinit),max(yinit)))
# plt.title("Initial positions of 2K particles")
# plt.show()

# #Tile the hexagonal shape
# lentile=max(xinit)-min(xinit)
# heightile=max(yinit)-min(yinit)
# xnumtiles=boxr/lentile
# ynumtiles=boxr/heightile
# print(lentile,heightile)

# lentile=max(xinit)-min(xinit) #Length of each tile
# heightile=max(yinit)-min(yinit) #height of each tile
# hexlayers=int(boxr/lentile) #Number of hexagonal layers

##Prepare modified B field
filehandler = open("modifiedfield.obj", 'rb')  ##2k particles 5.5hrs to run 2500 iterations just for settling down
gridr = pickle.load(filehandler)
gridz = pickle.load(filehandler)
Evalsr = pickle.load(filehandler)
Evalsz = pickle.load(filehandler)
rmax = pickle.load(filehandler)
separationsheath = pickle.load(filehandler)
separationhor1 = pickle.load(filehandler)
separationhor2 = pickle.load(filehandler)
firstpoint = pickle.load(filehandler)
filehandler.close()


def interpolate(r):
    x = r[0]
    y = r[1]
    r0 = numpy.sqrt(r[0] ** 2 + r[1] ** 2)
    z0 = r[2]
    if z0 > sheathd:
        return numpy.array([0, 0, 0])
    else:
        if r0 < rmax:
            Eresult = numpy.array([0, 0])
            indleft = (abs(r0) - abs(firstpoint)) / separationhor1
            indlow = (abs(z0) - abs(firstpoint)) / separationsheath
            s1 = isclose((indleft ** 10) ** (1.0 / 10), int(indleft))
            s2 = isclose((indlow ** 10) ** (1.0 / 10), int(indlow))
            # if s1==True and s2==True:
            # 	indlow=int(indlow)
            # 	indleft=int(indleft)
            # 	print("case 1")
            # 	gridEs=[[Evalsr[int(indlow)],Evalsz[int(indleft)]]]
            # 	gridpoints=[[gridr[int(indlow)],gridz[int(indleft)]]]
            # elif s1 ==True and s2 == False:
            # 	print("case 2")
            # 	indlow=int(indlow)
            # 	indleft=int(indleft)
            # 	gridEs=[[Evalsr[indlow][indleft],Evalsz[indlow][indleft]],[Evalsr[indlow+1][indleft],Evalsz[indlow+1][indleft]]]
            # 	gridpoints=[[gridr[indlow][indleft],gridz[indlow][indleft]],[gridr[indlow+1][indleft],gridz[indlow+1][indleft]]]

            # elif s1 ==False and s2 == True:
            # 	print("case 3")
            # 	indleft=int(indleft)
            # 	indlow=int(indlow)
            # 	gridEs=[[Evalsr[indlow][indleft],Evalsz[indlow][indleft]],[Evalsr[indlow][indleft+1],Evalsz[indlow][indleft+1]]]
            # 	gridpoints=[[gridr[indlow][indleft],gridz[indlow][indleft]],[gridr[indlow][indleft+1],gridz[indlow][indleft+1]]]

            # else:
            indlow = int(indlow)
            indleft = int(indleft)
            gridEs = [[Evalsr[indlow][indleft], Evalsz[indlow][indleft]],
                      [Evalsr[indlow + 1][indleft], Evalsz[indlow + 1][indleft]],
                      [Evalsr[indlow + 1][indleft + 1], Evalsz[indlow + 1][indleft + 1]],
                      [Evalsr[indlow][indleft + 1], Evalsz[indlow][indleft + 1]]]
            gridpoints = [[gridr[indlow][indleft], gridz[indlow][indleft]],
                          [gridr[indlow + 1][indleft], gridz[indlow + 1][indleft]],
                          [gridr[indlow + 1][indleft + 1], gridz[indlow + 1][indleft + 1]],
                          [gridr[indlow][indleft + 1], gridz[indlow][indleft + 1]]]

        elif r0 > rmax:
            Eresult = numpy.array([0, 0])
            indleft = (boxr - firstpoint - rmax) / separationhor2
            indlow = (z0 - firstpoint) / separationsheath
            s1 = isclose((indleft ** 3) ** (1.0 / 3), int(indleft))
            s2 = isclose((indlow ** 3) ** (1.0 / 3), int(indlow))
            indleft = int(indleft)
            indlow = int(indlow)
            if s1 != False or s2 != False:
                print("Case 4")
            gridEs = [[Evalsr[indlow][indleft], Evalsz[indlow][indleft]],
                      [Evalsr[indlow + 1][indleft], Evalsz[indlow + 1][indleft]],
                      [Evalsr[indlow + 1][indleft + 1], Evalsz[indlow + 1][indleft + 1]],
                      [Evalsr[indlow][indleft + 1], Evalsz[indlow][indleft + 1]]]
            gridpoints = [[gridr[indlow][indleft], gridz[indlow][indleft]],
                          [gridr[indlow + 1][indleft], gridz[indlow + 1][indleft]],
                          [gridr[indlow + 1][indleft + 1], gridz[indlow + 1][indleft + 1]],
                          [gridr[indlow][indleft + 1], gridz[indlow][indleft + 1]]]

        d1 = (abs(r0 - gridEs[0][0]))
        d4 = (abs(r0 - gridEs[0][1]))
        d2 = (abs(r0 - gridEs[2][0]))
        d3 = (abs(r0 - gridEs[2][1]))
        dhorsq = d1 ** 2 + d2 ** 2
        dvertsq = d3 ** 2 + d4 ** 2
        Etop = (d2 ** 2 / dhorsq) * numpy.array(gridEs[1]) + (d1 ** 2 / dhorsq) * numpy.array(gridEs[2])
        Ebottom = (d2 ** 2 / dhorsq) * numpy.array(gridEs[0]) + (d1 ** 2 / dhorsq) * numpy.array(gridEs[3])
        Efinal = (d3 ** 2 / dvertsq) * Ebottom + (d4 ** 2 / dvertsq) * Etop
        theta = numpy.arctan(abs(r[1] / r[0]))
        if r[0] > 0 and r[1] > 0:  # Quadrant 1
            pass
        elif r[0] < 0 and r[1] > 0:  # Quadrant 2
            theta = math.pi - theta
        elif r[0] < 0 and r[1] < 0:  # Quadrant 3
            theta += math.pi
        elif r[0] > 0 and r[1] < 0:  # Quadrant 4
            theta = 2 * math.pi - theta
        else:
            print("Problem: dust grain at x and y positions", [r[0], r[1]])
        Efinal = [Efinal[0] * numpy.cos(theta), Efinal[0] * numpy.sin(theta), Efinal[1]]
        return numpy.array([Efinal[0], Efinal[1], 0])


##Create dictionary of particles from pickle object
position = []
numparticles = 400
names = []
for i in numpy.arange(numparticles):
    names.append('g%s' % i)
dustdict = {name: Dust(md, radd, lambdaD, phia, Zd * e, [0, 0, 0], [0, 0, 0], [0, 0, 0]) for name in names}
for i, j in zip(dustdict, numpy.arange(len(dustdict))):
    dustdict[i].pos = initpositions[j]

# Create pairs
##Check every pair of particle interaction
list1 = list(dustdict.keys())
list2 = list1
pairlist = list(itertools.product(list1, list2))
pairs = set()
for x, y in pairlist:
    if x != y:
        pairs.add(frozenset((x, y)))
pairs = list(pairs)
for i in numpy.arange(len(pairs)):
    pairs[i] = list(pairs[i])

removelist = []
for i in pairs:
    if i[0] == i[1]:
        removelist.append(i)
pairs = [i for i in pairs if i not in removelist]

# Interact and iterate
iterationsB = 5000
inititerations = 100
g9velcheck = []
g9poscheck = []
g9acccheck = []

for i in tqdm(numpy.arange(inititerations)):
    pairsfinal = []
    for b in pairs:
        if dustdict[b[0]].intergraind(dustdict[b[1]]) == True:
            pairsfinal.append(b)
        else:
            pass  # pairsfinal.append(b)
    for j in pairsfinal:
        interactfield = dustdict[j[0]].selffield(dustdict[j[1]])
        dustdict[j[0]].selffieldmany(interactfield)
        dustdict[j[1]].selffieldmany(-interactfield)
    for k in dustdict:
        dustdict[k].steptest(numpy.array(dustdict[k].multifields))
        # acc.append(dustdict[k].getselfacc())
        # vel.append(dustdict[k].getselfvel())
        position.append(dustdict[k].getselfpos())

for l in dustdict:
    dustdict[l].Bswitch = True

for i in tqdm(numpy.arange(iterationsB)):
    pairsfinal = []
    for b in pairs:
        if dustdict[b[0]].intergraind(dustdict[b[1]]) == True:
            pairsfinal.append(b)
        else:
            pass  # pairsfinal.append(b)
    # g9velcheck.append(dustdict['g29'].getselfvel())
    # g9poscheck.append(dustdict['g29'].getselfpos())
    # g9acccheck.append(dustdict['g29'].getselfacc())
    ##Pick out the pairs that are less than 5 lambdaDs apart
    for j in pairsfinal:
        interactfield = dustdict[j[0]].selffield(dustdict[j[1]])
        dustdict[j[0]].selffieldmany(interactfield)
        dustdict[j[1]].selffieldmany(-interactfield)
    for k in dustdict:
        dustdict[k].steptest(numpy.array(dustdict[
                                             k].multifields))  # +interpolate(dustdict[k].getselfpos())) #Comment out interpolate function if want to turn of Gibson modified E field
        # acc.append(dustdict[k].getselfacc())
        # vel.append(dustdict[k].getselfvel())
        position.append(dustdict[k].getselfpos())

# Sort the positions of the particles
newx = [i[0] for i in position]
newy = [i[1] for i in position]
newz = [i[2] for i in position]

t = numpy.array(numpy.sort([list(numpy.arange(int(len(newx) / numparticles))) * numparticles]))[0]
df = pd.DataFrame({"time": t, "x": newx, "y": newy, "z": newz})


def update_graph(num):
    data = df[df['time'] == num]
    point.set_data(data.x, data.y)
    point.set_3d_properties(data.z)
    return point,


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=90., azim=90)
if min(newx) == max(newx):
    ax.set_xlim([-10 * lambdaD, 10 * lambdaD])
else:
    ax.set_xlim([min(newx), max(newx)])
if min(newy) == max(newy):
    ax.set_ylim([-10 * lambdaD, 10 * lambdaD])
else:
    ax.set_ylim([min(newy), max(newy)])
if min(newz) == max(newz):
    ax.set_zlim([-10 * lambdaD, 10 * lambdaD])
else:
    ax.set_zlim([min(newz) - 0.5 * min(newz), max(newz)])

ax.set_xlim([-rmax * 1.5, rmax * 1.5])
ax.set_ylim([-rmax * 1.5, rmax * 1.5])

data = df[df['time'] == 0]
point, = ax.plot(data.x, data.y, data.z, linestyle="", marker=".")
plt.cla()
# for i in [[0,0],[0,1],[1,0],[1,1]]:
# 	circx=((-1)**(i[0]))*numpy.arange(100)*boxr*0.9
# 	circy=((-1)**(i[1]))*numpy.sqrt(boxr**2-circx**2)
# 	ax.plot(circx,circy,'r-')
ani = matplotlib.animation.FuncAnimation(fig, update_graph, frames=iterationsB + inititerations, interval=1, blit=True)
#ani.save('RotationnoGibsvdrifttimes1.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
plt.show()
