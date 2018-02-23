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
from math import isclose

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
dt = 0.0000001
sheathd = 10 * lambdaD
electrodeV = abs((kb * Te / (2 * e)) * (numpy.log(2 * math.pi * me / mi)))  # potential at electrode
wallV = electrodeV  # cylindrical sides of wall same potential
radinfluence = 0.002  # 10*lambdaD
dipolea = boxr / 100.
mu0 = 4 * math.pi * 10 ** (-7)  # Permeaility free space
Bmom = ((2 * math.pi * (0.003) ** 3) * 0.014 / mu0) * numpy.array(
    [0, 0, 1])  # Nm/T #At 3mm away I want the B to be 0.014T
magBmom = numpy.sqrt(Bmom[0] ** 2 + Bmom[1] ** 2 + Bmom[2] ** 2)
Bmomhat = numpy.array(Bmom) / magBmom
dipolepos = [0, 0, -0.0005]

omegatau = 0.01  # Referenced by konopka in his experiments for ion collisions
omega = e * 0.014 / mi  # Tau_i (nu_i) is the one changing i.e. time between collisions (collisional frequency)


class ion:
    def __init__(self, pos, vel, acc):
        self.pos = numpy.array(pos)
        self.vel = numpy.array(vel)
        self.acc = numpy.array(acc)
        self.pos1 = numpy.array(pos)
        self.vel1 = numpy.array(vel)
        self.charge = e
        self.dr = 0

    def constB(self):
        return numpy.array([0, 0.014, 0])

    def constE(self, E=[0, 0, 1]):
        return numpy.array(E)

    def updateeulerforward(self, B, E=[0, 0, 1]):
        ##Euler forward - should be unstable!
        self.pos = self.pos + dt * self.vel
        self.vel = ((self.charge / mi) * numpy.cross(self.vel, B) + numpy.array(E) * (self.charge / mi)) * dt + self.vel

    def updateeulerback(self, B, E=[0, 0, 1]):  # Also unstable?
        self.pos = self.pos - dt * self.vel
        self.vel = self.vel - ((self.charge / mi) * numpy.cross(self.vel, B) * dt + numpy.array(E) * (self.charge / mi))

    def updateRK4(self, B, E=[0, 0, 1]):
        ##RK4 integration
        fv1 = (self.charge / mi) * numpy.cross(self.vel1, B) + numpy.array(E) * (self.charge / mi)
        fy1 = self.vel1
        v1 = self.vel1 + 0.5 * dt * fv1
        y1 = self.pos1 + 0.5 * dt * fy1
        fv2 = (self.charge / mi) * numpy.cross(v1, B) + numpy.array(E) * (self.charge / mi)
        fy2 = v1
        v2 = self.vel1 + 0.5 * dt * fv2
        y2 = self.pos1 + 0.5 * dt * fy2
        fv3 = (self.charge / mi) * numpy.cross(v2, B) + numpy.array(E) * (self.charge / mi)
        fy3 = v2
        v3 = self.vel1 + dt * fv3
        y3 = self.pos1 + dt * fy3
        fv4 = (self.charge / mi) * numpy.cross(v3, B) + numpy.array(E) * (self.charge / mi)
        fy4 = v3
        self.vel = self.vel1 + (1. / 6.) * dt * (fv1 + 2. * fv2 + 2. * fv3 + fv4)
        self.pos = self.pos1 + (1. / 6.) * dt * (fy1 + 2. * fy2 + 2. * fy3 + fy4)
        self.vel1 = self.vel
        self.dr = numpy.array(self.pos) - numpy.array(self.pos1)
        self.pos1 = self.pos

    def updateboris(self, B):
        pass

    def getselfvel(self):
        return self.vel

    def getselfpos(self):
        return self.pos


def runsim(iterations=2000):
    iterations = iterations
    vinit = numpy.sqrt(kb * Ti / mi)
    ion1 = ion(pos=[0, 0, 0], vel=[vinit, 0, 0], acc=[0, 0, 0])
    # ion2=ion(pos=[0,0,0],vel=[vinit,0,0],acc=[0,0,0])
    # ion3=ion(pos=[0,0,0],vel=[vinit,0,0],acc=[0,0,0])

    position = []
    velocity = []
    # velocity2=[]
    # position2=[]
    # position3=[]
    for i in numpy.arange(iterations):
        ion1.updateRK4(B=ion1.constB())
        # ion2.updateeulerforward(B=ion2.constB())
        # ion3.updateeulerback(B=ion3.constB())
        position.append(ion1.getselfpos())
        velocity.append(numpy.linalg.norm(ion1.getselfvel()))
    # velocity2.append(numpy.linalg.norm(ion2.getselfvel()))
    # position2.append(ion2.getselfpos())
    # position3.append(ion3.getselfpos())

    newx = [i[0] for i in position]
    newy = [i[1] for i in position]
    newz = [i[2] for i in position]
    # newx2=[i[0] for i in position2]
    # newy2=[i[1] for i in position2]
    # newz2=[i[2] for i in position2]
    # newx3=[i[0] for i in position3]
    # newy3=[i[1] for i in position3]
    # newz3=[i[2] for i in position3]
    direction = numpy.cross(ion1.constE(), ion1.constB()) / (
    numpy.linalg.norm(numpy.cross(ion1.constE(), ion1.constB())))
    return newx, newy, newz, numpy.array(direction), ion1


# newx,newy,newz,direct,ion1=runsim(iterations=20000)


def plotrunsim(newx, newy, newz):
    newx = newx
    newy = newy
    newz = newz
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(newx, newy, newz, 'r-', label="RK4")
    # ax.plot(newx2,newy2,newz2,'b-',label="Euler Forwards")
    # ax.plot(newx3,newy3,newz3,'m-',label="Euler backwards")

    plt.legend()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.show()


# (newx,newy,newz)

def checkdrift():
    iterations = numpy.arange(200, 20000, 500)
    particledrift = []
    for i in tqdm(iterations):
        x, y, z, direction, ion1 = runsim(i)
        EXBdrift = numpy.linalg.norm(ion1.constE()) / numpy.linalg.norm(ion1.constB())
        iterationtime = dt * i
        totald = EXBdrift * iterationtime
        theorydrift = numpy.linalg.norm(direction * totald)
        particlepos = numpy.linalg.norm(numpy.array([x[-1], y[-1], z[-1]]) * direction)
        particledrift.append([theorydrift, particlepos])
    particledrift = numpy.array(numpy.mat(particledrift).T)
    plt.figure()
    plt.scatter(iterations, particledrift[1] / particledrift[0])
    plt.ylabel("r_simulation/r_theory")
    plt.xlabel("Iteration number")
    plt.title("Ratio of simulated drift position vs theoretical drift position")
    plt.show()


# checkdrift()

def compareFEM(iterations=500):
    iterations = iterations
    vinit = numpy.sqrt(kb * Ti / mi)
    ion1 = ion(pos=[0, 0, 0], vel=[vinit, 0, 0], acc=[0, 0, 0])
    ion2 = ion(pos=[0, 0, 0], vel=[vinit, 0, 0], acc=[0, 0, 0])
    ion3 = ion(pos=[0, 0, 0], vel=[vinit, 0, 0], acc=[0, 0, 0])

    KERK4 = []
    KEEF = []
    KEEB = []
    for i in numpy.arange(iterations):
        ion1.updateRK4(B=ion1.constB(), E=[0, 0, 0])
        ion2.updateeulerforward(B=ion2.constB(), E=[0, 0, 0])
        ion3.updateeulerback(B=ion3.constB(), E=[0, 0, 0])
        KERK4.append(0.5 * mi * numpy.linalg.norm(ion1.getselfvel()) ** 2)
        KEEF.append(0.5 * mi * numpy.linalg.norm(ion2.getselfvel()) ** 2)
        KEEB.append(0.5 * mi * numpy.linalg.norm(ion3.getselfvel()) ** 2)

    # Compare Kinetic energy NB: when E field is [0,0,0] between RK4 an Euler forward
    plt.figure()
    plt.plot(numpy.arange(len(KERK4)), KERK4, label="RK4")
    plt.plot(numpy.arange(len(KEEB)), KEEB, 'r*', label="Euler backward")
    plt.plot(numpy.arange(len(KEEF)), KEEF, label="Euler forward")
    plt.legend()
    plt.ylim([0, max([max(KEEF) * 1.1, max(KEEB) * 1.1])])
    plt.title("No E field present, iterations = %s" % iterations)
    plt.xlabel("Iteration number")
    plt.ylabel("Kinetic energy of particle")
    plt.show()


# compareFEM()

def thermalkickexb(iterations, tau):
    iterations = iterations
    tau = tau
    tkick = int(tau / dt)  # Every tkick iterations give particle a thermal kick
    if (tau / dt) - tkick >= 0.5:
        tkick += 1  # Round upwards if decimal values >0.5
    if tkick == 0:
        raise ValueError("Error: dt>tau")
    vinit = numpy.sqrt(kb * Ti / mi)
    ion1 = ion(pos=[0, 0, 0], vel=[vinit, 0, 0], acc=[0, 0, 0])
    position = []
    velocity = []

    for i in tqdm(numpy.arange(int(iterations / tkick))):
        theta = math.pi * numpy.random.random_sample()
        phi = 2 * math.pi * numpy.random.random_sample()
        rhat = numpy.array([numpy.sin(theta) * numpy.cos(phi), numpy.sin(theta) * numpy.sin(phi), numpy.cos(theta)])
        ion1.vel1 = rhat * numpy.sqrt(kb * Ti / mi)
        ion1.updateRK4(B=ion1.constB())
        position.append(ion1.getselfpos())
        velocity.append(numpy.linalg.norm(ion1.getselfvel()))
        for j in numpy.arange(tkick - 1):
            ion1.updateRK4(B=ion1.constB())
            position.append(ion1.getselfpos())
            velocity.append(numpy.linalg.norm(ion1.getselfvel()))

    newx = [i[0] for i in position]
    newy = [i[1] for i in position]
    newz = [i[2] for i in position]
    return newx[-1] - newx[0]  ##Note need to return in whatever direction we're drifting!


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot(newx,newy,newz,'r--',label="RK4")
# ax.scatter(newx[-1],newy[-1],newz[-1],'ro',label="End",s=50)
# ax.scatter(newx[0],newy[0],newz[0],'bo',label="Start",s=50)
# plt.legend()
# ax.set_xlabel("x")
# ax.set_ylabel("y")
# ax.set_zlabel("z")
# plt.title("Thermal kick included at every dt=%s seconds"%tau)
# plt.show()

# drift = thermalkickexb(iterations=20000,tau=10**(-7))
# ion1=ion(pos=[0,0,0],vel=[numpy.sqrt(kb*Ti/mi),0,0],acc=[0,0,0])
# driftnocol=-dt*2*10**7*numpy.linalg.norm(ion1.constE())/numpy.linalg.norm(ion1.constB())

def averagekickeffect(iterations, tau, runs):
    driftdistancecol = []
    for i in tqdm(numpy.arange(runs), desc="run number for specific tau"):
        driftdistancecol.append(thermalkickexb(iterations, tau))
    ion1 = ion(pos=[0, 0, 0], vel=[numpy.sqrt(kb * Ti / mi), 0, 0], acc=[0, 0, 0])
    # driftnocol=-dt*iterations*numpy.linalg.norm(ion1.constE())/numpy.linalg.norm(ion1.constB())
    driftdistancecol = numpy.array(driftdistancecol)  # /driftnocol
    return driftdistancecol


# drifts=averagekickeffect(iterations = 10**7, tau = 10**(-5), runs = 1)
# fig = plt.figure()
# plt.scatter(numpy.arange(len(drifts)),drifts)
# plt.xlabel("Iteration number")
# plt.ylabel("driftlength_collision/driftlength_nocollision")
# plt.title("Plot of ratio of drift with thermal kicks to drift using pure EXB")
# fig.show()

def bootstrap(drifts, bsit=2000):  # Bootstrap iteration is to take bsit resamplings.
    drifts.sort()
    mean = numpy.mean(drifts)
    samplesind = numpy.random.choice(len(drifts), (bsit, len(drifts)))
    samples = []
    sampleavs = []
    for i in numpy.arange(bsit):
        sample = [drifts[k] for k in samplesind[i]]
        samples.append(sample)
        sampleavs.append(numpy.mean(sample))
    delta = [j - mean for j in sampleavs]
    delta.sort()
    tenpt = math.ceil(0.1 * len(drifts))
    nintypt = math.ceil(0.9 * (len(drifts)))
    error = [delta[tenpt - 1], delta[nintypt - 1]]
    return mean, error


# bs = bootstrap(drifts)

######################Compare different values of omega*tau and saving the drift values
# taulist = [0.01/omega, 0.05/omega, 0.1/omega, 0.5/omega]
# filehandler = open(b"drifts.obj",'wb')
# iterations=10**5
# runs = 2
# for i in tqdm(taulist,desc="Taus"):
# 	tau=i
# 	drifts=averagekickeffect(iterations = iterations, tau = tau, runs = runs)
# 	bs = bootstrap(drifts, bsit = 2000)
# 	pickle.dump(drifts,filehandler)
# 	pickle.dump(bs,filehandler)
# pickle.dump(numpy.array(taulist)*omega,filehandler)
# filehandler.close()

######################### Take a look at the saved object for drifts at different tau's
# filehandler = open("driftschangetautom.obj",'rb')
# drifts=[]
# bs=[]
# for i in numpy.arange(4):
# 	drifts.append(pickle.load(filehandler))
# 	bs.append(pickle.load(filehandler))
# omegataulist = pickle.load(filehandler)
# filehandler.close()

# So for omegatau is 0.01 the drift ratio after 1s is 0.0225916 (one run)
# So for omegatau is 0.01 the drift ratio after 2s is  0.01847764(one run)

# # ################ At different time steps for a fixed omega*tau
# iterationslist = [3/dt]
# filehandler = open(b"driftsconsttauwithEomegataupt275.obj",'wb')
# runs = 2
# for i in iterationslist:
# 	tau=0.275/(e*0.014/mi)
# 	drifts=averagekickeffect(iterations = i, tau = tau, runs = runs)
# 	bs = bootstrap(drifts,bsit = 2000)
# 	pickle.dump(drifts,filehandler)
# 	pickle.dump(bs,filehandler)
# pickle.dump(numpy.array(iterationslist),filehandler)
# filehandler.close()

# # ##############Take a look at the saved object for drifts with different times ####
# filehandler = open("driftsconsttauwithEomegataupt1.obj", 'rb')
# tau = 0.1/(e*0.14/mi)
# drifts = []
# drifts2 = []
# drifts3 = []
# bs = []
# bs2 =[]
# bs3 = []
#
# for i in numpy.arange(3):
#     drifts.append(pickle.load(filehandler))
#     bs.append(pickle.load(filehandler))
# iterationslist = pickle.load(filehandler)
# filehandler.close()
# timelist = numpy.array(iterationslist) * dt
# driftsav = [i[0] for i in bs]
# driftsav = numpy.array(driftsav) / (-dt * (1 / 0.014) * numpy.array(iterationslist))
# poserr = [i[1][0] for i in bs]
# negerr = [i[1][1] for i in bs]
# poserr = numpy.array(poserr) / (dt * (1 / 0.014) * numpy.array(iterationslist))
# negerr = numpy.array(negerr) / (dt * (1 / 0.014) * numpy.array(iterationslist))
#
# func = lambda n: (omega * tau) ** n / (1 + (omega * tau) ** n) - driftsav[-1]
# n = scipy.optimize.fsolve(func, 2)
#
# fig = plt.figure()
# plt.errorbar(timelist, driftsav, yerr=[poserr, negerr], fmt='o', ecolor='g', label='Average drifts')
# plt.xlabel("Run time(s)")
# plt.ylabel("Simualted drift/theoretical drift (m)")
# plt.legend()
# plt.title(
#     r"Different runtime, (omega*tau)^n for n=%s"%n)
# fig.show()

#Aggregate all the different drift ratios for the different tau simulations
filehandler = open("/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/driftsconsttauwithEomegataupt1.obj", 'rb')
filehandler2 = open("/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/driftsconsttauwithEtauE6.obj", 'rb')
filehandler3 = open("/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/driftsconsttauwithE.obj", 'rb')
filehandler4 = open("/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/driftsconsttauwithEomegataupt15.obj",'rb')
filehandler5 = open("/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/driftsconsttauwithEomegataupt25.obj",'rb')
filehandler5 = open("/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/driftsconsttauwithEomegataupt25.obj",'rb')
filehandler6 = open("/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/driftsconsttauwithEomegataupt3.obj",'rb')
filehandler7 = open("/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/driftsconsttauwithEomegataupt1real.obj",'rb')
filehandler8 = open("/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/driftsconsttauwithEomegataupt2.obj",'rb')
filehandler9 = open("/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/driftsconsttauwithEomegataupt07.obj",'rb')
filehandler10 = open("/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/driftsconsttauwithEomegataupt325.obj",'rb')
filehandler11 = open("/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/driftsconsttauwithEomegataupt275.obj",'rb')

# drifts = []
# drifts2 = []
# drifts3 = []
# drifts4= []
# drifts5 = []
# drifts6 = []
# drifts7= []
# drifts8 = []
# drifts9 = []
# drifts10 = []
# drifts11 = []
#
# bs = []
# bs2 =[]
# bs3 = []
# bs4 =[]
# bs5 =[]
# bs6 = []
# bs7 = []
# bs8 = []
# bs9 = []
# bs10 = []
# bs11 = []

driftobjects = [
    {'name': 'driftobject0', 'drifts': [], 'bs': [],  'num_of_runs': 1, 'filename': 'driftsconsttauwithEomegataupt275.obj'},
    {'name': 'driftobject1', 'drifts': [], 'bs': [],  'num_of_runs': 4, 'filename': 'driftsconsttauwithEtauE6.obj'},
    {'name': 'driftobject2', 'drifts': [], 'bs': [],  'num_of_runs': 2, 'filename': 'driftsconsttauwithE.obj'},
]

for driftobject in tqdm(driftobjects):
    filehandler = open(
        "/Users/yuewang/Dropbox/Msci-DustyPlasmas/Code/objects/{}".format(driftobject['filename']),
        'rb'
    )

    for i in numpy.arange(driftobject['num_of_runs']):
        driftobject['drifts'].append(pickle.load(filehandler))
        driftobject['bs'].append(pickle.load(filehandler))

    driftobject['iterations'] = pickle.load(filehandler)
    filehandler.close()
    driftobject['timelist'] = numpy.array(driftobject['iterations']) * dt
    driftobject['avx'] = [i[0] for i in driftobject['bs']]
    driftobject['av'] = numpy.array(driftobject['avx']) / (-dt * (1 / 0.014) * numpy.array(driftobject['timelist']) / dt)
    driftobject['poserr'] = [i[1][0] for i in driftobject['bs']]
    driftobject['negerr'] = [i[1][1] for i in driftobject['bs']]
    driftobject['poserr'] = numpy.array(driftobject['poserr']) / (dt * (1 / 0.014) * numpy.array(driftobject['timelist']) / dt)
    driftobject['negerr'] = numpy.array(driftobject['negerr']) / (dt * (1 / 0.014) * numpy.array(driftobject['timelist']) / dt)

# for i in numpy.arange(3):
#     drifts.append(pickle.load(filehandler))
#     bs.append(pickle.load(filehandler))
# iterationslist = pickle.load(filehandler)
# filehandler.close()
# timelist = numpy.array(iterationslist) * dt
# driftsavx = [i[0] for i in bs]
# driftsav = numpy.array(driftsavx) / (-dt * (1 / 0.014) * numpy.array([1/dt,2/dt,3/dt]))
# poserr = [i[1][0] for i in bs]
# negerr = [i[1][1] for i in bs]
# poserr = numpy.array(poserr) / (dt * (1 / 0.014) * numpy.array([1/dt,2/dt,3/dt]))
# negerr = numpy.array(negerr) / (dt * (1 / 0.014) * numpy.array([1/dt,2/dt,3/dt]))
#
# for i in numpy.arange(4):
#     drifts2.append(pickle.load(filehandler2))
#     bs2.append(pickle.load(filehandler2))
# iterationslist2 = pickle.load(filehandler2)
# filehandler2.close()
# timelist2 = numpy.array(iterationslist2) * dt/omega
# driftsav2x = [i[0] for i in bs2]
# driftsav2 = numpy.array(driftsav2x) / (-dt * (1 / 0.014) * numpy.array([0.5/dt,1/dt,2/dt,3/dt]))
# poserr2 = [i[1][0] for i in bs2]
# negerr2 = [i[1][1] for i in bs2]
# poserr2 = numpy.array(poserr2) / (dt * (1 / 0.014) * numpy.array([0.5/dt,1/dt,2/dt,3/dt]))
# negerr2 = numpy.array(negerr2) / (dt * (1 / 0.014) * numpy.array([0.5/dt,1/dt,2/dt,3/dt]))
#
# for i in numpy.arange(2):
#     drifts3.append(pickle.load(filehandler3))
#     bs3.append(pickle.load(filehandler3))
# iterationslist3 = pickle.load(filehandler3)
# filehandler3.close()
# timelist3 = numpy.array(iterationslist3) * dt/omega
# driftsav3x = [i[0] for i in bs3]
# driftsav3 = numpy.array(driftsav3x) / (-dt * (1 / 0.014) * numpy.array([0.5/dt,1/dt]))
# poserr3 = [i[1][0] for i in bs3]
# negerr3 = [i[1][1] for i in bs3]
# poserr3 = numpy.array(poserr3) / (dt * (1 / 0.014) * numpy.array([0.5/dt,1/dt]))
# negerr3 = numpy.array(negerr3) / (dt * (1 / 0.014) * numpy.array([0.5/dt,1/dt]))
#
# for i in numpy.arange(1):
#     drifts4.append(pickle.load(filehandler4))
#     bs4.append(pickle.load(filehandler4))
# iterationslist4 = pickle.load(filehandler4)
# filehandler4.close()
# timelist4 = numpy.array(iterationslist4) * dt
# driftsav4x = [i[0] for i in bs4]
# driftsav4 = numpy.array(driftsav4x) / (-dt * (1 / 0.014) * numpy.array([3/dt]))
# poserr4 = [i[1][0] for i in bs4]
# negerr4 = [i[1][1] for i in bs4]
# poserr4 = numpy.array(poserr4) / (dt * (1 / 0.014) * numpy.array([3/dt]))
# negerr4 = numpy.array(negerr4) / (dt * (1 / 0.014) * numpy.array([3/dt]))
#
# for i in numpy.arange(1):
#     drifts5.append(pickle.load(filehandler5))
#     bs5.append(pickle.load(filehandler5))
# iterationslist5 = pickle.load(filehandler5)
# filehandler5.close()
# timelist5 = numpy.array(iterationslist4) * dt
# driftsav5x = [i[0] for i in bs5]
# driftsav5 = numpy.array(driftsav5x) / (-dt * (1 / 0.014) * numpy.array([3/dt]))
# poserr5 = [i[1][0] for i in bs5]
# negerr5 = [i[1][1] for i in bs5]
# poserr5 = numpy.array(poserr5) / (dt * (1 / 0.014) * numpy.array([3/dt]))
# negerr5 = numpy.array(negerr5) / (dt * (1 / 0.014) * numpy.array([3/dt]))
#
# for i in numpy.arange(1):
#     drifts6.append(pickle.load(filehandler6))
#     bs6.append(pickle.load(filehandler6))
# iterationslist6 = pickle.load(filehandler6)
# filehandler6.close()
# timelist6 = numpy.array(iterationslist6) * dt
# driftsav6x = [i[0] for i in bs6]
# driftsav6 = numpy.array(driftsav6x) / (-dt * (1 / 0.014) * numpy.array([3/dt]))
# poserr6 = [i[1][0] for i in bs6]
# negerr6 = [i[1][1] for i in bs6]
# poserr6 = numpy.array(poserr6) / (dt * (1 / 0.014) * numpy.array([3/dt]))
# negerr6 = numpy.array(negerr6) / (dt * (1 / 0.014) * numpy.array([3/dt]))
#
# for i in numpy.arange(1):
#     drifts7.append(pickle.load(filehandler7))
#     bs7.append(pickle.load(filehandler7))
# iterationslist7 = pickle.load(filehandler7)
# filehandler7.close()
# timelist7 = numpy.array(iterationslist7) * dt
# driftsav7x = [i[0] for i in bs7]
# driftsav7 = numpy.array(driftsav7x) / (-dt * (1 / 0.014) * numpy.array([3/dt]))
# poserr7 = [i[1][0] for i in bs7]
# negerr7 = [i[1][1] for i in bs7]
# poserr7 = numpy.array(poserr7) / (dt * (1 / 0.014) * numpy.array([3/dt]))
# negerr7 = numpy.array(negerr7) / (dt * (1 / 0.014) * numpy.array([3/dt]))
#
#
# for i in numpy.arange(1):
#     drifts8.append(pickle.load(filehandler8))
#     bs8.append(pickle.load(filehandler8))
# iterationslist8 = pickle.load(filehandler8)
# filehandler8.close()
# timelist8 = numpy.array(iterationslist8) * dt
# driftsav8x = [i[0] for i in bs8]
# driftsav8 = numpy.array(driftsav8x) / (-dt * (1 / 0.014) * numpy.array([3/dt]))
# poserr8 = [i[1][0] for i in bs8]
# negerr8 = [i[1][1] for i in bs8]
# poserr8 = numpy.array(poserr8) / (dt * (1 / 0.014) * numpy.array([3/dt]))
# negerr8 = numpy.array(negerr8) / (dt * (1 / 0.014) * numpy.array([3/dt]))
#
#
# for i in numpy.arange(1):
#     drifts9.append(pickle.load(filehandler9))
#     bs9.append(pickle.load(filehandler9))
# iterationslist9 = pickle.load(filehandler9)
# filehandler9.close()
# timelist9 = numpy.array(iterationslist9) * dt
# driftsav9x = [i[0] for i in bs9]
# driftsav9 = numpy.array(driftsav9x) / (-dt * (1 / 0.014) * numpy.array([3/dt]))
# poserr9 = [i[1][0] for i in bs9]
# negerr9 = [i[1][1] for i in bs9]
# poserr9 = numpy.array(poserr9) / (dt * (1 / 0.014) * numpy.array([3/dt]))
# negerr9 = numpy.array(negerr9) / (dt * (1 / 0.014) * numpy.array([3/dt]))
#
#
# for i in numpy.arange(1):
#     drifts10.append(pickle.load(filehandler10))
#     bs10.append(pickle.load(filehandler10))
# iterationslist10 = pickle.load(filehandler10)
# filehandler10.close()
# timelist10 = numpy.array(iterationslist10) * dt
# driftsav10x = [i[0] for i in bs10]
# driftsav10 = numpy.array(driftsav10x) / (-dt * (1 / 0.014) * numpy.array([3/dt]))
# poserr10 = [i[1][0] for i in bs10]
# negerr10 = [i[1][1] for i in bs10]
# poserr10 = numpy.array(poserr10) / (dt * (1 / 0.014) * numpy.array([3/dt]))
# negerr10 = numpy.array(negerr10) / (dt * (1 / 0.014) * numpy.array([3/dt]))
#
# for i in numpy.arange(1):
#     drifts11.append(pickle.load(filehandler11))
#     bs11.append(pickle.load(filehandler11))
# iterationslist11 = pickle.load(filehandler11)
# filehandler11.close()
# timelist11 = numpy.array(iterationslist10) * dt
# driftsav11x = [i[0] for i in bs11]
# driftsav11 = numpy.array(driftsav11x) / (-dt * (1 / 0.014) * numpy.array([3/dt]))
# poserr11 = [i[1][0] for i in bs11]
# negerr11 = [i[1][1] for i in bs11]
# poserr11 = numpy.array(poserr11) / (dt * (1 / 0.014) * numpy.array([3/dt]))
# negerr11 = numpy.array(negerr11) / (dt * (1 / 0.014) * numpy.array([3/dt]))
#
# import matplotlib
# font = {'family' : 'normal',
#         'weight' : 'bold',
#         'size'   : 25}
#
# matplotlib.rc('font', **font)
# fig=plt.figure()
# finalratios = [driftsav[-1],driftsav2[-1],driftsav3[-1],driftsav4[-1],driftsav5[-1],driftsav6[-1],driftsav7[-1],driftsav8[-1],driftsav9[-1],driftsav10[-1],driftsav11[-1]]
# omegatau = [0.01,0.033813824472513944,0.3381382447251395, 0.15,0.25,0.3,0.1,0.2,0.07,0.325,0.275]
# plt.errorbar(omegatau,finalratios,yerr=[poserr3[-1]*numpy.array(omegatau)*1.5/(numpy.linalg.norm(omegatau)),negerr3[-1]*1.5*numpy.array(omegatau)/numpy.linalg.norm(omegatau)] ,fmt='o',capthick=2)
# theoryx=numpy.arange(0,max(omegatau),10**(-2))
# plt.plot(theoryx, theoryx**3/(1+theoryx)**3,'r-',label = r'$y = \frac{(\omega*\tau)^3}{1+(\omega*\tau)^3}$')
# # plt.plot(theoryx, theoryx**2/(1+theoryx)**2, label='n=2')
# plt.xlabel(r"$\omega*\tau$",fontsize = 25)
# plt.ylabel("simulated drift velocity / theoretical drift velocity", fontsize = 25)
# plt.title("Reduction factor for ion-drift velocity \n due to ion-neutral collisions")
# plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
# plt.legend(loc=2)
# from pylab import rcParams
# rcParams['figure.figsize'] = 16, 10
# # plt.show()
# plt.savefig('fig', format='eps')