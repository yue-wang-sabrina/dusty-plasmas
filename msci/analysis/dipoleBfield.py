import numpy
import matplotlib.pyplot as plt
import msci.analysis.constants as const
import math
from scipy import special
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3

from IPython import get_ipython

ipython = get_ipython()
ipython.magic('load_ext autoreload')
ipython.magic('autoreload 2')


def pointdipoleB(r):
    r = numpy.array(r) - numpy.array(const.dipolepos)
    magr = numpy.sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2)
    rhat = numpy.array(r) / magr
    magBmom = numpy.sqrt(const.Bmom[0] ** 2 + const.Bmom[1] ** 2 + const.Bmom[2] ** 2)
    Bmomhat = numpy.array(const.Bmom) / magBmom
    return (const.mu0 * (3 * numpy.dot(rhat, numpy.array(const.Bmom)) * rhat - const.Bmom) / (4 * math.pi * magr ** 3))


def finitesolenoid(r, N=50, L=0.5, rad=0.1,
                   I=1000):  # N coils over length of solenoid L, rad is radius of solenoid cross section, current I
    magr = numpy.sqrt(r[0] ** 2 + r[1] ** 2)
    sigmaplus = r[2] + L / 2.
    sigmaminus = r[2] - L / 2.
    n = N / L
    ksquarep = 4 * const.radd * magr / (sigmaplus ** 2 + (rad + magr) ** 2)
    kp = numpy.sqrt(ksquarep)
    ksquarem = 4 * const.radd * magr / (sigmaminus ** 2 + (rad + magr) ** 2)
    km = numpy.sqrt(ksquarem)
    phip = numpy.arctan(abs(sigmaplus / (rad - magr)))
    phim = numpy.arctan(abs(sigmaminus / (rad - magr)))

    def zeta(phi, ksquare):
        return special.ellipeinc(phi, ksquare) - special.ellipe(ksquare) * special.ellipkinc(phi,
                                                                                             ksquare) / special.ellipk(
            ksquare)

    def heuman(phi, ksquare):
        return (special.ellipkinc(phi, (1. - ksquare)) / special.ellipk((1 - ksquare))) + (
                                                                                          2. / math.pi) * special.ellipk(
            ksquare) * zeta(phi, (1 - ksquare))

    Br = (const.mu0 * n * I / math.pi) * numpy.sqrt(const.radd / magr) * (
        (((2 - ksquarep) / (2 * kp)) * special.ellipk(kp) - special.ellipe(kp) / kp)
        - (((2 - ksquarem) / (2 * km)) * special.ellipk(km) - special.ellipe(km) / km))
    x, y = ((numpy.array(r) / magr) * Br)[0:2]
    Bz = (const.mu0 * n * I / 4.) * ((sigmaplus * kp * special.ellipk(kp) / (math.pi * numpy.sqrt(magr * rad)) + (
        rad - magr) * sigmaplus * heuman(phip, kp) / abs((rad - magr) * sigmaplus)) - (
                                         sigmaminus * km * special.ellipk(km) / (math.pi * numpy.sqrt(magr * rad)) + (
                                             rad - magr) * sigmaminus * heuman(phim, km) / abs(
                                             (rad - magr) * sigmaminus)))
    return numpy.array([x, y, 0])

def plotBfield(): #3D plot of B field
    heightbox = const.sheathd
    widthbox = 0.002
    dx = 0.
    stepsize = widthbox / 3
    stepsizeheight = heightbox / 3

    # X, Y, Z = numpy.meshgrid(numpy.arange(-0.5, 0.5 , stepsize),numpy.arange(-0.5, 0.5 , stepsize), numpy.arange(-0.5, 0.5, stepsize))
    X, Y, Z = numpy.meshgrid(numpy.arange(-widthbox, widthbox, stepsize), numpy.arange(-widthbox, widthbox, stepsize),
                             numpy.arange(-heightbox + dx, heightbox + dx, stepsizeheight))

    U = []
    V = []
    W = []
    for i in numpy.arange(len(Y)):
        rowU = []
        rowV = []
        rowW = []
        for j in numpy.arange(len(Y[i])):
            rowU2 = []
            rowV2 = []
            rowW2 = []
            for k in numpy.arange(len(Y[i][j])):
                if abs(numpy.sqrt(X[i][j][k] ** 2 + Y[i][j][k] ** 2 + Z[i][j][k] ** 2)) < 0.05 * const.boxz:
                    rowU2.append(0)
                    rowV2.append(0)
                    rowW2.append(0)
                else:
                    newpoint = pointdipoleB([X[i][j][k], Y[i][j][k], Z[i][j][k]])
                    rowU2.append(newpoint[0])
                    rowV2.append(newpoint[1])
                    rowW2.append(newpoint[2])
            rowU.append(rowU2)
            rowV.append(rowV2)
            rowW.append(rowW2)
        U.append(numpy.array(rowU))
        V.append(numpy.array(rowV))
        W.append(numpy.array(rowW))
    U = numpy.array(U)
    V = numpy.array(V)
    W = numpy.array(W)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    # ax.view_init(elev=90., azim=90)
    Ulist = abs(numpy.hstack(numpy.hstack(U)))
    Vlist = abs(numpy.hstack(numpy.hstack(V)))
    Wlist = abs(numpy.hstack(numpy.hstack(W)))
    maximumarrow = max(max(Ulist), max(Vlist), max(Wlist))
    Q = ax.quiver3D(X, Y, Z, U, V, W, length=0.0001, arrow_length_ratio=0.5,
                    normalize=True)  # length is just a multiplication factor
    # P = ax.scatter(X,Y,Z)
    # maxx=max(Ulist)
    # maxy=max(Vlist)
    # maxz=max(Wlist)
    # ax.set_xlim([-widthbox,widthbox])
    # ax.set_ylim([-widthbox,widthbox])
    # ax.set_zlim([-heightbox,heightbox+dx])
    plt.title("3D dipole Magnetic field")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    fig.show()


plotBfield()