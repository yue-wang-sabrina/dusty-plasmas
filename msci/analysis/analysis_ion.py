import numpy
import math
import matplotlib.pyplot as plt
import pickle

from tqdm import tqdm

import msci.analysis.constants as const
from msci.particles.ions import Ion
from msci.utils.utils import bootstrap

dt = 0.0000001


class IonAnalysis:
    def __init__(self, tau, omega):
        self.position = None
        self.iterations = None
        self.directiondrift = None
        self.method = None
        self.ion = None
        self.particledrift = None
        self.const = const
        self.KERK4 = None
        self.KEEF = None
        self.KEEB = None
        self.drift = None
        self.driftdistancecol = None
        self.tau = tau
        self.omega = omega
        self.dt = Ion(pos=[0, 0, 0], vel=[numpy.sqrt(const.kb * const.Ti / const.mi), 0, 0], acc=[0, 0, 0],
                      omega=self.omega, tau=self.tau).dt

    def runsim(self, iterations, method):  # Run a single ion motion through E and B field and make a plot
        self.iterations = iterations
        self.method = method
        vinit = numpy.sqrt(const.kb * const.Ti / const.mi)
        self.ion = Ion(pos=[0, 0, 0], vel=[vinit, 0, 0], acc=[0, 0, 0], omega=self.omega, tau=self.tau)

        position = []
        velocity = []

        if method == 'RK4':
            for i in numpy.arange(self.iterations):
                self.ion.updateRK4(B=self.ion.constB())
                position.append(self.ion.getselfpos())
                velocity.append(numpy.linalg.norm(self.ion.getselfvel()))

        elif method == 'EF':
            for i in numpy.arange(self.iterations):
                self.ion.updateeulerforward(B=self.ion.constB())
                position.append(self.ion.getselfpos())
                velocity.append(numpy.linalg.norm(self.ion.getselfvel()))

        elif method == 'EB':
            for i in numpy.arange(self.iterations):
                self.ion.updateeulerback(B=self.ion.constB())
                position.append(self.ion.getselfpos())
                velocity.append(numpy.linalg.norm(self.ion.getselfvel()))

        elif method == 'BORIS':
            for i in numpy.arange(self.iterations):
                self.ion.updateboris(B=self.ion.constB())
                position.append(self.ion.getselfpos())
                velocity.append(numpy.linalg.norm(self.ion.getselfvel()))

        else:
            raise ValueError("{} does not exist".format(method))

        self.position = numpy.array(position)
        self.directiondrift = numpy.cross(self.ion.constE(), self.ion.constB()) / (
            numpy.linalg.norm(numpy.cross(self.ion.constE(), self.ion.constB())))

    def plotrunsim(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(self.position[:, 0], self.position[:, 1], self.position[:, 2], 'r-', label="{}".format(self.method))
        plt.legend()
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.set_zlabel("z (m)")
        fig.show()

    def checkdrift(self, method, iterations=numpy.arange(200, 1000, 500)):
        drift = []
        for i in tqdm(iterations):
            self.runsim(i, method=method)
            EXBdrift = numpy.linalg.norm(self.ion.constE()) / numpy.linalg.norm(self.ion.constB())
            iterationtime = dt * i
            totald = EXBdrift * iterationtime
            theorydrift = numpy.linalg.norm(self.directiondrift * totald)
            particlepos = numpy.linalg.norm(
                numpy.array([self.position[-1, 0], self.position[-1, 1], self.position[-1, 2]]) * self.directiondrift
            )
            drift.append([theorydrift, particlepos])
        self.iterations = iterations
        self.particledrift = numpy.array(numpy.mat(drift).T)

    def plotcheckdrift(self):
        plt.figure()
        plt.scatter(self.iterations, self.particledrift[1] / self.particledrift[0])
        plt.ylabel("r_simulation/r_theory")
        plt.xlabel("Number of iterations")
        plt.title("Ratio of simulated drift position vs theoretical drift position")
        plt.show()

    def compareFEM(self, iterations):
        self.iterations = iterations
        vinit = numpy.sqrt(const.kb * const.Ti / const.mi)
        ion1 = Ion(pos=[0, 0, 0], vel=[vinit, 0, 0], acc=[0, 0, 0], omega=self.omega, tau=self.tau)
        ion2 = Ion(pos=[0, 0, 0], vel=[vinit, 0, 0], acc=[0, 0, 0], omega=self.omega, tau=self.tau)
        ion3 = Ion(pos=[0, 0, 0], vel=[vinit, 0, 0], acc=[0, 0, 0], omega=self.omega, tau=self.tau)

        KERK4 = []
        KEEF = []
        KEEB = []

        for i in numpy.arange(iterations):
            ion1.updateRK4(B=ion1.constB(), E=[0, 0, 0])
            ion2.updateeulerforward(B=ion2.constB(), E=[0, 0, 0])
            ion3.updateeulerback(B=ion3.constB(), E=[0, 0, 0])
            KERK4.append(0.5 * const.mi * numpy.linalg.norm(ion1.getselfvel()) ** 2)
            KEEF.append(0.5 * const.mi * numpy.linalg.norm(ion2.getselfvel()) ** 2)
            KEEB.append(0.5 * const.mi * numpy.linalg.norm(ion3.getselfvel()) ** 2)

        self.KERK4 = KERK4
        self.KEEF = KEEF
        self.KEEB = KEEB

    def plotcompareFEM(self):
        # Compare Kinetic energy NB: when E field is [0,0,0] between RK4 an Euler forward
        plt.figure()
        plt.plot(numpy.arange(len(self.KERK4)), self.KERK4, 'b-', label="RK4")
        plt.plot(numpy.arange(len(self.KEEB)), self.KEEB, 'r-', label="Euler backward")
        plt.plot(numpy.arange(len(self.KEEF)), self.KEEF, 'g-', label="Euler forward")
        plt.legend()
        plt.ylim([0, max(self.KERK4) * 1000])
        plt.title("No E field present, iterations = %s" % self.iterations)
        plt.xlabel("Iteration number")
        plt.ylabel("Kinetic energy of particle")
        plt.show()

    def thermalkickexb(self, iterations):
        iterations = iterations
        tau = self.tau
        tkick = int(tau / dt)  # Every tkick iterations give particle a thermal kick
        if (tau / dt) - tkick >= 0.5:
            tkick += 1  # Round upwards if decimal values >0.5
        if tkick == 0:
            raise ValueError("Error: dt>tau")
        vinit = numpy.sqrt(const.kb * const.Ti / const.mi)
        ion1 = Ion(pos=[0, 0, 0], vel=[vinit, 0, 0], acc=[0, 0, 0], omega=self.omega, tau=self.tau)
        position = []
        velocity = []

        for i in tqdm(numpy.arange(int(iterations / tkick))):
            theta = math.pi * numpy.random.random_sample()
            phi = 2 * math.pi * numpy.random.random_sample()
            rhat = numpy.array([numpy.sin(theta) * numpy.cos(phi), numpy.sin(theta) * numpy.sin(phi), numpy.cos(theta)])
            ion1.vel1 = rhat * numpy.sqrt(const.kb * const.Ti / const.mi)
            ion1.updateRK4(B=ion1.constB())
            position.append(list(ion1.getselfpos()))
            velocity.append(numpy.linalg.norm(ion1.getselfvel()))
            for j in numpy.arange(tkick - 1):
                ion1.updateRK4(B=ion1.constB())
                position.append(list(ion1.getselfpos()))
                velocity.append(numpy.linalg.norm(ion1.getselfvel()))
        direction = numpy.cross(ion1.constE(), ion1.constB())
        drift = numpy.dot(numpy.array(position[-1]) - numpy.array(position[0]), direction)
        self.drift = drift
        self.position = position

    def averagekickeffect(self, iterations, runs):
        tau = self.tau
        driftdistancecol = []
        for p in tqdm(numpy.arange(runs), desc="Doing another run for tau={}".format(tau)):
            self.thermalkickexb(iterations)
            driftdistancecol.append(self.drift)
        ion1 = Ion(pos=[0, 0, 0], vel=[numpy.sqrt(const.kb * const.Ti / const.mi), 0, 0], acc=[0, 0, 0],
                   omega=self.omega, tau=self.tau)
        # driftnocol=-dt*iterations*numpy.linalg.norm(ion1.constE())/numpy.linalg.norm(ion1.constB())
        self.driftdistancecol = numpy.array(driftdistancecol)  # /driftnocol

    def changetau(self, taulist, iterations, runs, filename, security=False):
        if security:
            filehandler = open("{}".format(filename), 'wb')
            for i in tqdm(taulist, desc="Taus"):
                tau = i
                self.averagekickeffect(iterations=iterations, runs=runs)
                bs = bootstrap(self.driftdistancecol, bsit=2000)
                pickle.dump(self.driftdistancecol, filehandler)
                pickle.dump(bs, filehandler)
            pickle.dump(numpy.array(taulist) * self.omega, filehandler)  # Dumps the values of omega*tau into pickle
            pickle.dump(iterations, filehandler)
            filehandler.close()
        else:
            print("Security is False, make it true to change the file {}".format(filename))

    def changetime(self, tau, iterationlist, runs, filename, security=False):
        if security:
            print("Altering file for drifts")
            filehandler = open("{}".format(filename), 'wb')
            for i in iterationlist:
                tau = tau
                self.averagekickeffect(iterations=i, runs=runs)
                bs = bootstrap(self.driftdistancecol, bsit=2000)
                pickle.dump(self.driftdistancecol, filehandler)
                pickle.dump(bs, filehandler)
            pickle.dump(numpy.array(iterationlist) * self.dt, filehandler)
            filehandler.close()
        else:
            print("Security is False, make it true to change the file {}".format(filename))
