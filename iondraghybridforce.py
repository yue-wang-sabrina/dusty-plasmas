##This file is mainly for looking at B field effects of crystals that are already formed.
import numpy
import scipy
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

# Some constants
radd = 1.5 * 10 ** (-6)  # radius of dust particle
e = 1.60217662 * 10 ** (-19)  # electron charge
e0 = 8.85418782 * 10 ** (-12)  # perimittivity free space
Te = 46000  # Electrons usually a few EV #Non-thermal plasma, E field couples with electron more efficiently, higher temperature
Ti = 310  # Ion thermal temperature, room temperature
kb = 1.38064852 * 10 ** (-23)  # boltzmans constant
ne0 = 1. * 10 ** (15)  # initial electron density
ni0 = ne0  # Quasi-neutral initial conditions
md = 1000. * (4 / 3) * math.pi * radd ** 3  # Mass of dust particle
mi = 39.948 * 1.66053904 * 10 ** (-27)  # argon mass ion
me = 9.10938356 * 10 ** (-31)  # electron mass
cs = numpy.sqrt(kb * Te / mi)  # Bohm speed = sound speed
lambdade = (kb * Te * e0 / (ne0 * e ** 2)) ** 0.5  # Electron debye rad
lambdadi = (kb * Ti * e0 / (ni0 * e ** 2)) ** 0.5  # Ion debye rad
lambdaD = 1. / (1. / (lambdadi ** 2) + 1. / (lambdade ** 2)) ** 0.5  # dusty debye radius
boxr = 523 * lambdaD  # cylinder radius
boxz = 0.001  # cylinder height
g = -9.8  # gravitational acceleration
dt = 0.0001  # Much longer than charging time which is order nanoseconds
sheathd = 10 * lambdaD
electrodeV = abs((kb * Te / (2 * e)) * (numpy.log(2 * math.pi * me / mi)))  # potential at electrode
wallV = electrodeV  # cylindrical sides of wall same potential
radinfluence = 20 * lambdaD
dipolea = boxr / 100.
mu0 = 4 * math.pi * 10 ** (-7)  # Permeaility free space
Bmom = ((2 * math.pi * (0.003) ** 3) * 0.014 / mu0) * numpy.array(
    [0, 0, 1])  # Nm/T #At 1cm away I want the B to be 0.014T
magBmom = numpy.sqrt(Bmom[0] ** 2 + Bmom[1] ** 2 + Bmom[2] ** 2)
Bmomhat = numpy.array(Bmom) / magBmom
dipolepos = [0, 0, -0.0005]


def force(u, Zd=-9388.3579633332938):
    vT = numpy.sqrt(kb * Ti / mi)
    mach = numpy.array([u, 0, 0])
    machmag = numpy.linalg.norm(mach)
    LAMBDA = numpy.sqrt(1. / (numpy.exp(-machmag ** 2 / 2.) * lambdadi ** (-2) + lambdade ** (-2)))
    beta = abs(Zd * e ** 2 / (LAMBDA * Ti * kb))
    z = abs(Zd) * e ** 2 / (4 * math.pi * e0 * radd * Te * kb)
    tau = Te / Ti
    coloumblog = 5.
    force = numpy.sqrt(2 * math.pi) * radd ** 2 * ni0 * mi * vT ** 2 * \
            (numpy.sqrt(math.pi / 2) * special.erf(u / numpy.sqrt(2)) * \
             (1 + u ** 2 + (1 - u ** (-2)) * (1 + 2 * z * tau) + 4 * z ** 2 * tau ** 2 * u ** (-2) * numpy.log(
                 coloumblog)) + \
             (u ** (-1) * (1 + 2 * z * tau + u ** 2 - 4 * z ** 2 * tau ** 2 * numpy.log(coloumblog)) * numpy.exp(
                 -u ** 2 / 2.)))
    return force


ulist = numpy.arange(0.1, 100, 0.1)
acceleration = []
for i in ulist:
    acceleration.append(numpy.linalg.norm(force(i)))

plt.loglog(ulist, acceleration)
plt.xlabel("u = ion flow velocity / ion thermal velocity")
plt.ylabel("Ion drag force")
plt.title("Ion drag force vs normalised ion flow velocity")
plt.show()
