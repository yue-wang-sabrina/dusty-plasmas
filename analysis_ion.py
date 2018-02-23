import numpy
import math
import matplotlib.pyplot as plt

from tqdm import tqdm

import constants as const
from particles.ions import Ion

dt = 0.0000001


class IonAnalysis:
    def __init__(self):
        self.position = None
        self.iterations = None
        self.direction = None
        self.method = None
        self.ion = None
        self.particledrift = None
        self.const = const
        self.KERK4 = None
        self.KEEF = None
        self.KEEB = None

    def runsim(self, iterations, method):  # Run a single ion motion through E and B field and make a plot
        self.iterations = iterations
        self.method = method
        vinit = numpy.sqrt(const.kb * const.Ti / const.mi)
        self.ion = Ion(pos=[0, 0, 0], vel=[vinit, 0, 0], acc=[0, 0, 0])

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
        self.direction = numpy.cross(self.ion.constE(), self.ion.constB()) / (
            numpy.linalg.norm(numpy.cross(self.ion.constE(), self.ion.constB())))

    def plotrunsim(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(self.position[:, 0], self.position[:, 1], self.position, 'r-', label="{}".format(self.method))
        plt.legend()
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.set_zlabel("z (m)")
        fig.show()

    def checkdrift(self, method, iterations=numpy.arange(200, 20000, 500)):
        drift = []
        for i in tqdm(iterations):
            self.runsim(i, method=method)
            EXBdrift = numpy.linalg.norm(self.ion.constE()) / numpy.linalg.norm(self.ion.constB())
            iterationtime = dt * i
            totald = EXBdrift * iterationtime
            theorydrift = numpy.linalg.norm(self.direction * totald)
            particlepos = numpy.linalg.norm(
                numpy.array([self.position[-1, 0], self.position[-1, 1], self.position[-1, 2]]) * self.direction
            )
            drift.append([theorydrift, particlepos])
        self.particledrift = numpy.array(numpy.mat(drift).T)

    def plotcheckdrift(self):
        plt.figure()
        plt.scatter(self.iterations, self.particledrift[1] / self.particledrift[0])
        plt.ylabel("r_simulation/r_theory")
        plt.xlabel("Iteration number")
        plt.title("Ratio of simulated drift position vs theoretical drift position")
        plt.show()

    def compareFEM(self, iterations):
        self.iterations = iterations
        vinit = numpy.sqrt(const.kb * const.Ti / const.mi)
        ion1 = Ion(pos=[0, 0, 0], vel=[vinit, 0, 0], acc=[0, 0, 0])
        ion2 = Ion(pos=[0, 0, 0], vel=[vinit, 0, 0], acc=[0, 0, 0])
        ion3 = Ion(pos=[0, 0, 0], vel=[vinit, 0, 0], acc=[0, 0, 0])

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
        plt.plot(numpy.arange(len(self.KERK4)), self.KERK4, label="RK4")
        plt.plot(numpy.arange(len(self.KEEB)), self.KEEB, 'r*', label="Euler backward")
        plt.plot(numpy.arange(len(self.KEEF)), self.KEEF, label="Euler forward")
        plt.legend()
        plt.ylim([0, max([max(self.KEEF) * 1.1, max(self.KEEB) * 1.1])])
        plt.title("No E field present, iterations = %s" % self.iterations)
        plt.xlabel("Iteration number")
        plt.ylabel("Kinetic energy of particle")
        plt.show()

